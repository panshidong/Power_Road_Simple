"""Physical three-phase OpenDSS operation on TAMU's complete regional circuits.

Two operators are available (``power.operator``):

* ``uniform_grid`` (legacy, 2026-10-06 runs): one load multiplier for a whole
  region, lowered along a fixed grid until every element is within its
  nameplate rating. Because the TAMU base case already loads some substation
  transformers and 69 kV lines at 2-4.7x nameplate, this sheds 60-80 % of a
  healthy region and makes repairs able to *reduce* served load.
* ``local``: the normal rating of an element is max(nameplate, healthy
  base-case loading), so pre-existing synthetic-data overloads are not blamed
  on damage and the healthy grid serves its nominal load. Damaged states are
  checked against emergency ratings, ``emergency_factor`` x normal (TAMU's
  own transformer convention is EmergHKVA = 1.5 x kVA); a derated damaged
  transformer gets only ``factor`` x normal, without emergency loading. Any
  violation is removed by shedding only the loads downstream of the violated
  element. Tie switches may pick up dark load only as far as it fits without
  shedding customers who were already supplied.

Both are feasible heuristics, not an OPF or guaranteed optimum.
"""
from __future__ import annotations
import concurrent.futures as futures
import math, multiprocessing, os, signal, subprocess, sys, tempfile, time
import numpy as np
from collections import OrderedDict, defaultdict
from pathlib import Path
from .common import RUNTIME, atomic_json, cached_json, digest, local_path, read_json, rows
from .progress import phase, report
from .native_env import child_environment


class PowerInfeasible(RuntimeError):
    pass


class NativeACFailure(RuntimeError):
    def __init__(self, message, *, returncode, directory, retryable):
        super().__init__(message)
        self.returncode=returncode;self.directory=Path(directory);self.retryable=retryable


def physical_differences(actual, expected, path=""):
    """Strict physical identity, with only the previously validated float tolerance.

    Timing/counter/diagnostic metadata is excluded. Signals, topology, trial
    outcomes, limits, served loads and every other physical field are compared.
    """
    if isinstance(expected,dict):
        if not isinstance(actual,dict):return [path+": type mismatch"]
        keys=set(expected)-{"performance"}
        if set(actual)-{"performance"}!=keys:return [path+": keys differ"]
        return [item for key in sorted(keys) for item in physical_differences(actual[key],expected[key],path+"/"+key)]
    if isinstance(expected,list):
        if not isinstance(actual,list) or len(actual)!=len(expected):return [path+": list length/type differs"]
        return [item for i,(a,b) in enumerate(zip(actual,expected)) for item in physical_differences(a,b,path+f"/{i}")]
    if isinstance(expected,float):
        if (not isinstance(actual,(float,int)) or isinstance(actual,bool) or
            not math.isfinite(actual) or not math.isfinite(expected) or
            not math.isclose(actual,expected,rel_tol=1e-9,abs_tol=1e-8)):
            return [path+": numeric mismatch"]
    elif type(actual) is not type(expected) or actual!=expected:return [path+": value/type mismatch"]
    return []


def downstream_loads(starts, incidence, pterm, term_first, terminals, loads_at_bus, skip, eps=1e-3):
    """Loads supplied through ``starts`` along the solved active-power direction.

    ``incidence`` maps a bus to (element, terminal) pairs, ``pterm`` holds the
    real power entering each element terminal (kW, positive into the element).
    From a bus, an element is followed only if power enters it there, and only
    to its terminals where power leaves. In a radial feeder this is the subtree
    below the violated element; in a mesh it is every load the flow can reach.
    """
    seen=set(starts); stack=list(starts); found=[]
    while stack:
        bus=stack.pop(); found.extend(loads_at_bus.get(bus,()))
        for element,terminal in incidence.get(bus,()):
            if element==skip or pterm[term_first[element]+terminal]<=eps:continue
            for other,next_bus in enumerate(terminals[element]):
                if other!=terminal and pterm[term_first[element]+other]< -eps and next_bus not in seen:
                    seen.add(next_bus); stack.append(next_bus)
    return found


_regional_engine = None


def _init_regional_worker(cfg, catalog, fingerprint, scratch, verbose):
    global _regional_engine
    if hasattr(os,"setsid"): os.setsid()
    for name in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","NUMEXPR_NUM_THREADS"):
        os.environ[name]="1"
    _regional_engine=PowerEngine(cfg,catalog,fingerprint,Path(scratch)/f"region-worker-{os.getpid()}",verbose=verbose)


def _solve_regional_state(region, damage, fresh=False):
    if fresh:
        _regional_engine._loaded_region=None
        return _regional_engine._compute(region,damage)
    return _regional_engine._regional_state(region,damage)


class PowerEngine:
    def __init__(self, cfg, catalog, fingerprint, scratch, *, region_workers=1, verbose=False):
        from opendssdirect import dss
        self.cfg=cfg["power"]; self.catalog=catalog; self.fingerprint=fingerprint
        self._full_cfg=cfg
        self.verbose=verbose; self.region_workers=max(1,min(region_workers,len(catalog["regions"])))
        self._pool=None
        self._pool_args=(cfg,catalog,fingerprint,str(scratch),verbose)
        self.prepared=local_path(cfg,"prepared"); self.cache=local_path(cfg,"cache")/fingerprint/"power"
        self.scratch=Path(scratch); self.scratch.mkdir(parents=True,exist_ok=True)
        self.dss=dss.NewContext()
        self.version=self.dss.Basic.Version()
        self._loaded_region=None; self._baseline=None
        self._changed_enabled={}; self._changed_terminals={}
        self.memory=OrderedDict()
        self._changed_loads=set(); self._static={}; self._references={}
        self.stats=defaultdict(float)
        self.loads=defaultdict(dict)
        for r in rows(self.prepared/"load_catalog.csv.gz"):
            r["kw"]=float(r["kw"]); self.loads[r["region"]][r["load_id"]]=r
        self.signal_loads=defaultdict(list)
        for r in catalog["signals"]:
            region=catalog["substation_region"][r["substation_id"]]
            self.signal_loads[(region,r["load_id"].lower())].append(r["signal_id"])

    def _regional_state(self, region, damage):
        key=digest(dict(region=region,damage=damage,engine=self.version,settings=self.cfg))
        if key in self.memory:
            self.stats["memory_hits"]+=1
            self.memory.move_to_end(key)
            return self.memory[key]
        path=self.cache/region/(key+".json")
        if path.exists():self.stats["disk_hits"]+=1
        with phase(f"AC region {region}; damaged_assets={len(damage)}",enabled=self.verbose):
            result=cached_json(path,lambda:self._compute(region,damage))
        self.memory[key]=result
        if len(self.memory)>128:self.memory.popitem(last=False)
        return result

    def fresh_region(self, region, damage):
        # Use an existing regional process, without creating another large native
        # circuit in the coordinator. This deliberately bypasses the state cache.
        if self._pool is None:
            self._loaded_region=None
            return self._compute(region,damage)
        try:return self._pool.submit(_solve_regional_state,region,damage,True).result()
        except BaseException:
            self.close(cancel=True)
            raise

    def close(self, *, cancel=False):
        pool,self._pool=self._pool,None
        if hasattr(self,"dss"):
            self.dss.Basic.ClearAll()
            self._loaded_region=None; self._baseline=None
        if pool is None:return
        if cancel:
            for child in list((getattr(pool,"_processes",None) or {}).values()):
                try:os.killpg(child.pid,signal.SIGTERM)
                except ProcessLookupError:child.terminate()
                except PermissionError:child.terminate()
        pool.shutdown(wait=True,cancel_futures=cancel)

    def _compute(self, region, damage):
        if not getattr(self,"cfg",{}).get("isolate_states",False):return self._operate(region,damage)
        started=time.monotonic()
        try:return self._operate_isolated(region,damage)
        except NativeACFailure as failure:
            if not failure.retryable or self.cfg.get("native_recovery_attempts",0)!=1:raise
            # Exactly one recovery, then an independent confirmation. Neither
            # run uses cache; any second fault, timeout, or mismatch propagates.
            audit=failure.directory/"recovery.json"
            record=dict(status="running",region=region,damage=damage,initial_returncode=failure.returncode,
                        initial_failure=str(failure),maximum_recovery_attempts=1,independent_confirmation_required=True)
            atomic_json(audit,record)
            report(f"AC {region}: native exit {failure.returncode}; one fresh recovery plus independent physical confirmation; {audit}")
            try:
                recovered=self._operate_isolated(region,damage)
                record["recovery_directory"]=str(self._last_isolated_directory);atomic_json(audit,record)
                confirmed=self._operate_isolated(region,damage)
                record["confirmation_directory"]=str(self._last_isolated_directory)
                differences=physical_differences(confirmed,recovered)
                if differences:
                    record["differences"]=differences
                    raise RuntimeError(f"{region}: independent AC recovery results disagree; no cache write; {audit}")
                self.stats["native_recovered_states"]+=1
                elapsed=time.monotonic()-started
                record.update(status="recovered_and_confirmed",physical_differences=[],wall_seconds=elapsed)
                atomic_json(audit,record)
                recovered["performance"].update(native_recovered_states=1,recovery_wall_seconds=elapsed,recovery_audit=str(audit))
                self._mark("recovered_and_confirmed",region=region,damage=damage,audit=str(audit))
                return recovered
            except BaseException as exc:
                record.update(status="failed",error=f"{type(exc).__name__}: {exc}",wall_seconds=time.monotonic()-started)
                atomic_json(audit,record)
                raise

    def _operate_isolated(self, region, damage):
        """A cache miss owns one fresh native process, including all trial scales.

        Native circuits never survive across regional states. The caller retains
        the normal cache lock; only a clean exit with a matching feasible result
        can populate that cache. Recovery, if enabled, is handled by _compute;
        this single attempt never retries or substitutes physics.
        """
        timeout=float(self.cfg["state_timeout_seconds"])
        if not math.isfinite(timeout) or timeout<=0:raise ValueError("Invalid AC state timeout")
        folder=Path(tempfile.mkdtemp(prefix=f"ac-{region}-",dir=self.scratch)).resolve()
        self._last_isolated_directory=folder
        request=dict(config=self._full_cfg,fingerprint=self.fingerprint,region=region,damage=damage,
                     directory=str(folder),verbose=self.verbose)
        token=digest(request);request_path=folder/"request.json";atomic_json(request_path,request)
        env=dict(os.environ,PYTHONPATH=str(RUNTIME))
        for name in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","NUMEXPR_NUM_THREADS"):
            env[name]="1"
        env=child_environment(self._full_cfg,env)
        started=time.monotonic();code=None;reason=None
        self.stats["isolated_calls"]+=1
        with (folder/"native.log").open("w") as log:
            # Inherit this worker's process group so supervisor cancellation
            # still owns the native child. This child creates no further solvers.
            with subprocess.Popen([sys.executable,"-u","-m","austin_runtime.ac_worker",str(request_path)],
                stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,cwd=RUNTIME,env=env) as child:
                self._mark("isolated_running",region=region,damage=damage,child_pid=child.pid,directory=str(folder))
                try:code=child.wait(timeout=timeout)
                except subprocess.TimeoutExpired:
                    child.kill();code=child.wait();reason=f"AC state exceeded {timeout:g}s hard timeout"
                except BaseException:
                    child.kill();child.wait();raise
        elapsed=time.monotonic()-started;self.stats["isolated_wall_seconds"]+=elapsed
        try:
            if reason:raise RuntimeError(reason)
            if code:raise RuntimeError(f"Native AC process exited {code}")
            envelope=read_json(folder/"result.json");result=envelope["result"]
            if (envelope["request_digest"]!=token or envelope["fingerprint"]!=self.fingerprint or
                envelope.get("native_allocator_sha256")!=(self.cfg.get("native_allocator_sha256") or None) or
                result["region"]!=region or result["damage"]!=damage or result["engine"]!=self.version or
                not all(result.get(k) is True for k in ("feasible","converged","controls_settled"))):
                raise ValueError("Native AC result identity or physical acceptance mismatch")
            for key,value in result["performance"].items():self.stats[key]+=value
            result["performance"].update(process_wall_seconds=elapsed,child_peak_rss_mib=envelope["peak_rss_mib"])
            self._mark("complete",region=region,damage=damage,child_pid=child.pid,directory=str(folder),
                       performance=result["performance"])
            atomic_json(folder/"status.json",dict(status="passed",returncode=code,wall_seconds=elapsed))
            return result
        except Exception as exc:
            diagnostic=dict(status="failed",returncode=code,wall_seconds=elapsed,region=region,damage=damage,
                            error=f"{type(exc).__name__}: {exc}",request_digest=token)
            marker=folder/"dss/ac_current.json"
            if marker.exists():diagnostic["last_ac_operation"]=read_json(marker)
            atomic_json(folder/"failure.json",diagnostic)
            self._mark("isolated_failed",**diagnostic,directory=str(folder))
            retryable=reason is None and code in (-signal.SIGILL,-signal.SIGABRT,-signal.SIGBUS,-signal.SIGFPE,-signal.SIGSEGV)
            raise NativeACFailure(f"{region}: {exc}; AC diagnostics retained in {folder}",
                                  returncode=code,directory=folder,retryable=retryable) from exc

    def evaluate(self, remaining):
        tasks=[]
        for region in self.catalog["regions"]:
            damage={key:factor for key,factor in remaining.items() if key.startswith("power:") and self.catalog["assets"][key]["region"]==region}
            tasks.append((region,damage))
        if self.region_workers==1:
            results=[self._regional_state(region,damage) for region,damage in tasks]
        else:
            if self._pool is None:
                report(f"AC regional pool: {self.region_workers} processes for {len(tasks)} complete original circuits")
                self._pool=futures.ProcessPoolExecutor(max_workers=self.region_workers,
                    mp_context=multiprocessing.get_context("spawn"),initializer=_init_regional_worker,initargs=self._pool_args)
            completed={}
            try:
                pending={self._pool.submit(_solve_regional_state,region,damage):region for region,damage in tasks}
                for future in futures.as_completed(pending):
                    region=pending[future];completed[region]=future.result()
                    if self.verbose:report(f"AC regions complete: {len(completed)}/{len(tasks)}; region={region}")
            except BaseException:
                self.close(cancel=True)
                raise
            # Completion order must not change summation order or output identities.
            results=[completed[region] for region,_ in tasks]
        zones=defaultdict(float); subs=defaultdict(float); signals={}
        for r in results:
            for k,v in r["zone_served_kw"].items():zones[k]+=v
            for k,v in r["substation_served_kw"].items():subs[k]+=v
            signals.update(r["signal_powered"])
        return dict(served_kw=sum(r["served_kw"] for r in results),nominal_kw=sum(r["nominal_kw"] for r in results),
            zone_served_kw=dict(zones),substation_served_kw=dict(subs),signal_powered=signals,
            regions=[{k:v for k,v in r.items() if k not in ("zone_served_kw","substation_served_kw","signal_powered")} for r in results])

    def _capture_baseline(self):
        """Capture mutable state of the frozen TAMU snapshot model, before faults.

        Unsupported element classes take the original cold-compile path. This is
        intentionally not a general OpenDSS state serializer (e.g. storage).
        """
        d=self.dss
        supported={"vsource","line","transformer","load","regcontrol","capacitor","capcontrol","fuse"}
        if any(name.split(".",1)[0].lower() not in supported for name in d.Circuit.AllElementNames()):
            return None
        taps={}
        for r in d.RegControls:
            name=r.Transformer().lower().removeprefix("transformer."); winding=r.TapWinding()
            d.Transformers.Name(name); d.Transformers.Wdg(winding)
            taps[(name,winding)]=d.Transformers.Tap()
        return dict(taps=taps,caps=[(c.Name(),list(c.States())) for c in d.Capacitors],
                    fuses=[(f.Name(),list(f.State())) for f in d.Fuses])

    def _mark(self, stage, **detail):
        # Survives a native abort/SIGKILL; heartbeats alone cannot identify the
        # native operation which was active when a worker disappeared.
        atomic_json(self.scratch/"ac_current.json",dict(stage=stage,pid=os.getpid(),unix_time=time.time(),**detail))

    def _restore_baseline(self):
        d=self.dss; baseline=self._baseline
        # Enable precisely the objects changed by fault isolation, then restore
        # controller internals AND actuator states. Reset alone does not undo taps.
        for name,enabled in self._changed_enabled.items():
            d.Circuit.SetActiveElement(name); d.CktElement.Enabled(enabled)
        self._changed_enabled.clear()
        d("Reset")
        d.CtrlQueue.ClearQueue(); d.CtrlQueue.ClearActions()
        for (name,winding),tap in baseline["taps"].items():
            d.Transformers.Name(name); d.Transformers.Wdg(winding); d.Transformers.Tap(tap)
        for name,states in baseline["caps"]:
            d.Capacitors.Name(name); d.Capacitors.States(states)
        for name,states in baseline["fuses"]:
            d.Fuses.Name(name); d.Fuses.State(states)
        for name,states in self._changed_terminals.items():
            d.Circuit.SetActiveElement(name)
            for terminal,conductor,opened in states:
                (d.CktElement.Open if opened else d.CktElement.Close)(terminal,conductor)
        self._changed_terminals.clear()
        # Local shedding edits individual loads; return them to their nominal P/Q.
        for name in getattr(self,"_changed_loads",()):
            meta=self.loads[self._loaded_region][name]
            d.Loads.Name(name); d.Loads.kW(meta["kw"]); d.Loads.kvar(float(meta["kvar"]))
        self._changed_loads=set()
        # Recreate the original healthy zero-load voltage initialization before
        # applying faults. Do not start from the preceding solved voltages/taps.
        d.Solution.LoadMult(1.)
        d("CalcVoltageBases")

    def _load_circuit(self, region, scale, closed):
        d=self.dss
        d.Basic.AllowChangeDir(False)
        d.Basic.DataPath(str(self.scratch))
        master=self.prepared/"models"/region/"Master_runtime.dss"
        reuse=(getattr(self,"_loaded_region",None)==region and self._baseline is not None
               and self.cfg.get("reuse_circuit",True))
        started=time.monotonic()
        self._mark("reset" if reuse else "load",region=region,scale=scale,closed=closed)
        if reuse:
            try:self._restore_baseline()
            except BaseException:
                self._loaded_region=None
                raise
            self.stats["circuit_resets"]+=1
            self.stats["reset_seconds"]+=time.monotonic()-started
            return
        self._loaded_region=None
        with phase(f"AC {region}: load circuit / voltage bases; load_scale={scale}; ties={len(closed)}",
                   enabled=getattr(self,"verbose",False)):
            d(f'Redirect "{master}"')
        self._changed_enabled={}; self._changed_terminals={}; self._changed_loads=set()
        self._baseline=self._capture_baseline() if getattr(self,"cfg",{}).get("reuse_circuit",True) else None
        self._loaded_region=region
        if hasattr(self,"stats"):
            self.stats["circuit_loads"]+=1
            self.stats["load_seconds"]+=time.monotonic()-started

    def _disable(self, name):
        d=self.dss
        d.Circuit.SetActiveElement(name)
        if d.CktElement.Name().lower()!=name.lower():raise ValueError("Component not found: "+name)
        self._changed_enabled.setdefault(name,d.CktElement.Enabled())
        d("Disable "+name)

    def _compile(self, region, damage, scale, closed, multipliers=None):
        self._load_circuit(region,scale,closed)
        # Element ratings/topology are read from the restored healthy circuit, before isolation.
        if getattr(self,"cfg",{}).get("operator","uniform_grid")=="local":self._region_static(region)
        d=self.dss
        d.Basic.DataPath(str(self.scratch))
        d("Set mode=snapshot controlmode=static maxcontroliter=100 maxiterations=100")
        d.Solution.LoadMult(scale)
        disabled=set(); deratings={}
        for key,factor in damage.items():
            for name in self.catalog["assets"][key]["isolation_transformers"]:
                if factor==0:
                    self._disable(name); disabled.add(name.split(".",1)[1].lower())
                else: deratings[name]=factor
        # Only disable a regulator controlling an isolated transformer. Cross-feeder
        # references elsewhere remain intact; no control object is deleted.
        for name in d.RegControls.AllNames():
            if name.lower()=="none":continue
            d.RegControls.Name(name)
            if d.RegControls.Transformer().lower().removeprefix("transformer.") in disabled: self._disable("RegControl."+name)
        for name in closed:
            d.Circuit.SetActiveElement(name)
            if d.CktElement.Name().lower()!=name.lower():raise ValueError("Switch not found: "+name)
            if not d.CktElement.Enabled():raise ValueError("Cannot close disabled component: "+name)
            self._changed_terminals[name]=[(t,c,d.CktElement.IsOpen(t,c)) for t in (1,2)
                                          for c in range(1,d.CktElement.NumConductors()+1)]
            for terminal in (1,2): d.CktElement.Close(terminal,0)
        expected=self.loads[region]
        if d.Loads.Count()!=len(expected):
            raise ValueError(f"{region}: native load count {d.Loads.Count()} differs from catalog {len(expected)}")
        for name,m in (multipliers or {}).items():
            meta=expected[name]
            d.Loads.Name(name); d.Loads.kW(meta["kw"]*m); d.Loads.kvar(float(meta["kvar"])*m)
            self._changed_loads.add(name)
        return deratings,disabled

    def _measure(self, region, scale, closed, deratings):
        d=self.dss; p=self.cfg
        voltage=dict(zip((n.lower() for n in d.Circuit.AllNodeNames()),map(float,d.Circuit.AllBusMagPu())))
        if not voltage or not all(math.isfinite(v) and v>=0 for v in voltage.values()):raise PowerInfeasible("Nonfinite voltages")
        zones=defaultdict(float); subs=defaultdict(float); powered={}
        actual_total=served=nominal=0.0; dark=0; energized_values=[]; missing=set(self.loads[region])
        for load in d.Loads:
            name=load.Name().lower(); meta=self.loads[region].get(name)
            if meta is None:raise ValueError(f"Unexpected regional load {region}/{name}")
            missing.remove(name); kw=meta["kw"]; nominal+=kw
            if not math.isclose(load.kW(),kw,rel_tol=1e-8,abs_tol=1e-7) or not math.isclose(load.kvar(),float(meta["kvar"]),rel_tol=1e-8,abs_tol=1e-7):
                raise ValueError(f"Native nominal P/Q differs from frozen load catalog: {region}/{name}")
            bus=d.CktElement.BusNames()[0].split(".")[0].lower()
            phases=[int(n) for n in d.CktElement.NodeOrder() if n]
            volts=[voltage.get(f"{bus}.{n}",0.0) for n in phases]
            energized_values.extend(v for v in volts if v>1e-4)
            energized=bool(volts) and all(v>1e-4 for v in volts)
            power=d.CktElement.Powers(); nc=d.CktElement.NumConductors()
            actual=sum(float(x) for x in power[:2*nc:2])
            if not math.isfinite(actual):raise PowerInfeasible("Nonfinite load power")
            actual_total+=actual
            delivered=max(0.0,min(kw,actual)) if energized else 0.0
            served+=delivered; zones[meta["zone"]]+=delivered; subs[meta["substation_id"]]+=delivered
            dark+=not energized
            good=energized and min(volts)>=p["voltage_min_pu"]-p["voltage_tolerance"] and max(volts)<=p["voltage_max_pu"]+p["voltage_tolerance"]
            good=good and delivered/max(kw,1e-12)>=p["signal_min_served_fraction"]
            for sid in self.signal_loads.get((region,name),[]):powered[sid]=bool(good)
        if missing:raise ValueError("Native load identities missing from solve")
        vmin=min(energized_values) if energized_values else None
        vmax=max(energized_values) if energized_values else None
        max_loading=0.0; unrated=[]
        for line in d.Lines:
            if not d.CktElement.Enabled():continue
            amps=[float(x) for x in d.CktElement.CurrentsMagAng()[::2]]
            current=max(amps,default=0.0)
            rating=float(d.Properties.Value("normamps"))
            if not math.isfinite(current):raise PowerInfeasible("Nonfinite line currents")
            if rating<=0:
                if current>1e-6:unrated.append(d.CktElement.Name())
            else:max_loading=max(max_loading,current/rating)
        for transformer in d.Transformers:
            if not d.CktElement.Enabled():continue
            name=d.CktElement.Name().lower(); values=list(map(float,d.CktElement.Powers()))
            nc=d.CktElement.NumConductors()
            for winding in range(1,transformer.NumWindings()+1):
                transformer.Wdg(winding); rating=transformer.kVA()*deratings.get(name,1.0)
                start=(winding-1)*2*nc; powers=values[start:start+2*nc]
                kva=math.hypot(sum(powers[::2]),sum(powers[1::2]))
                if not math.isfinite(kva):raise PowerInfeasible("Nonfinite transformer power")
                if rating<=0:
                    if kva>1e-6:unrated.append(name)
                else:max_loading=max(max_loading,kva/rating)
        source_kw=-float(d.Circuit.TotalPower()[0]); losses_kw=float(d.Circuit.Losses()[0])/1000
        residual=source_kw-actual_total-losses_kw
        balance_ok=abs(residual)<=max(.1,abs(source_kw)*p["balance_relative_tolerance"])
        feasible=(not unrated and balance_ok and max_loading<=1+p["thermal_tolerance"] and
            (vmin is None or vmin>=p["voltage_min_pu"]-p["voltage_tolerance"]) and
            (vmax is None or vmax<=p["voltage_max_pu"]+p["voltage_tolerance"]))
        return dict(region=region,feasible=bool(feasible),converged=True,controls_settled=True,
            load_scale=scale,served_kw=served,nominal_kw=nominal,actual_load_kw=actual_total,source_kw=source_kw,
            losses_kw=losses_kw,active_balance_residual_kw=residual,min_loaded_voltage_pu=vmin,max_loaded_voltage_pu=vmax,
            max_thermal_loading_ratio=max_loading,unrated_energized_elements=unrated,dark_loads=dark,
            zone_served_kw=dict(zones),substation_served_kw=dict(subs),signal_powered=powered,closed_switches=list(closed),
            voltage_limits=[p["voltage_min_pu"],p["voltage_max_pu"]],engine=self.version)

    def _feasible(self, region, damage, closed):
        attempts=[]
        for scale in self.cfg["load_scales"]:
            deratings,disabled=self._compile(region,damage,scale,closed)
            try:
                started=time.monotonic()
                self._mark("solve",region=region,scale=scale,damage=damage,closed=closed)
                with phase(f"AC {region}: power flow; load_scale={scale}; ties={len(closed)}",enabled=self.verbose):
                    self.dss.Solution.Solve()
                self.stats["solve_seconds"]+=time.monotonic()-started
                self.stats["solve_calls"]+=1
                if not self.dss.Solution.Converged() or not self.dss.Solution.ControlActionsDone():
                    raise PowerInfeasible("Power flow or controller iteration did not converge")
                for name in disabled:
                    self.dss.Circuit.SetActiveElement("Transformer."+name)
                    if self.dss.CktElement.Name().lower()!="transformer."+name or self.dss.CktElement.Enabled():
                        raise ValueError("Fault isolation failed: "+name)
                started=time.monotonic()
                self._mark("measure",region=region,scale=scale,damage=damage,closed=closed)
                with phase(f"AC {region}: check voltages / thermal limits / balance",enabled=self.verbose):
                    result=self._measure(region,scale,closed,deratings)
                self.stats["measure_seconds"]+=time.monotonic()-started
                if self.verbose:report(f"AC {region}: load_scale={scale}, feasible={result['feasible']}")
                result["isolated_transformers"]=sorted(disabled)
                result["rating_derating_factors"]=deratings
                attempts.append({k:result[k] for k in ("load_scale","feasible","min_loaded_voltage_pu","max_loaded_voltage_pu","max_thermal_loading_ratio","active_balance_residual_kw","unrated_energized_elements")})
                if result["feasible"]:
                    result["curtailment_trials"]=attempts
                    return result
            except self.dss.DSSException as exc: attempts.append(dict(load_scale=scale,error=str(exc)))
            except PowerInfeasible as exc: attempts.append(dict(load_scale=scale,error=str(exc)))
        raise PowerInfeasible(f"{region}: no feasible operating point on curtailment grid; trials={attempts}")

    def _tie_candidates(self, already_closed):
        d=self.dss; candidates=[]
        # Collect line properties before activating buses (active-object APIs are stateful).
        switches=[]
        for line in d.Lines:
            name=d.CktElement.Name()
            if name in already_closed or not d.CktElement.Enabled():continue
            if d.Properties.Value("switch").strip().lower() not in ("yes","true","1"):continue
            if not any(d.CktElement.IsOpen(t,0) for t in (1,2)):continue
            switches.append((name,[b.split(".")[0] for b in d.CktElement.BusNames()]))
        for name,buses in switches:
            live=[]
            for bus in buses[:2]:
                d.Circuit.SetActiveBus(bus)
                values=d.Bus.puVmagAngle()[::2]
                live.append(bool(values) and max(values)>.1)
            if len(live)==2 and live[0]!=live[1]:candidates.append(name)
        return sorted(candidates)

    def _operate(self, region, damage):
        if self.cfg.get("operator","uniform_grid")=="local":return self._operate_local(region,damage)
        started=time.monotonic(); before=dict(self.stats)
        best=self._feasible(region,damage,[]); evaluations=[]
        for iteration in range(self.cfg["max_switch_closures"]):
            # Reestablish the best state before identifying energized/dark endpoints.
            if iteration:best=self._feasible(region,damage,best["closed_switches"])
            candidates=self._tie_candidates(best["closed_switches"])[:self.cfg["max_switch_trials"]]
            incumbent=best
            for name in candidates:
                try:
                    trial=self._feasible(region,damage,best["closed_switches"]+[name])
                    evaluations.append(dict(switch=name,served_kw=trial["served_kw"],feasible=True))
                    if trial["served_kw"]>incumbent["served_kw"]+1e-6:incumbent=trial
                except PowerInfeasible as exc:evaluations.append(dict(switch=name,feasible=False,error=str(exc)))
            if incumbent is best:break
            best=incumbent
        best["switch_trials"]=evaluations
        best["damage"]=damage
        best["operator"]="descending uniform-curtailment grid plus limited greedy existing-tie closure; no optimality claim"
        best["performance"]={k:v-before.get(k,0) for k,v in self.stats.items()}
        best["performance"]["wall_seconds"]=time.monotonic()-started
        self._mark("complete",region=region,damage=damage,performance=best["performance"])
        return best

    # ------------------------------------------------------------------
    # Local operator: healthy-base-case limits and downstream shedding.
    # ------------------------------------------------------------------

    def _region_static(self, region):
        """PD element names, terminal buses and nameplate ratings, plus load buses.

        Read from the loaded circuit before any fault is applied, once per
        process and region. Names keep the native order of the batch arrays.
        """
        d=self.dss
        names=[n.lower() for n in d.PDElements.AllNames()]
        cached=self._static.get(region)
        if cached is not None and cached["names"]==names:return cached
        nt=np.asarray(d.PDElements.AllNumTerminals(),dtype=np.int64)
        nc=np.asarray(d.PDElements.AllNumConductors(),dtype=np.int64)
        if len(nt)!=len(names) or len(nc)!=len(names) or np.any(nt*nc<1):raise ValueError(f"{region}: unexpected PD element layout")
        width=nt*nc
        elem_start=np.concatenate(([0],np.cumsum(width)[:-1]))
        term_first=np.concatenate(([0],np.cumsum(nt)[:-1]))
        term_elem=np.repeat(np.arange(len(names)),nt)
        term_start=np.repeat(elem_start,nt)+(np.arange(int(nt.sum()))-np.repeat(term_first,nt))*np.repeat(nc,nt)
        index={n:i for i,n in enumerate(names)}
        kind=np.zeros(len(names),dtype=np.int8);terminals=[]
        for i,name in enumerate(names):
            d.Circuit.SetActiveElement(name)
            buses=[b.split(".")[0].lower() for b in d.CktElement.BusNames()][:int(nt[i])]
            terminals.append(buses+[buses[-1]]*(int(nt[i])-len(buses)))
            kind[i]=1 if name.startswith("line.") else 2 if name.startswith("transformer.") else 0
        line_rating=np.full(len(names),np.nan)
        for _ in d.Lines:line_rating[index[d.CktElement.Name().lower()]]=float(d.Lines.NormAmps())
        term_rating=np.full(int(nt.sum()),np.nan)
        for transformer in d.Transformers:
            i=index[d.CktElement.Name().lower()]
            for winding in range(1,min(transformer.NumWindings(),int(nt[i]))+1):
                transformer.Wdg(winding);term_rating[term_first[i]+winding-1]=float(transformer.kVA())
        incidence=defaultdict(list)
        for i,buses in enumerate(terminals):
            for t,bus in enumerate(buses):incidence[bus].append((i,t))
        load_nodes={};loads_at_bus=defaultdict(list)
        for load in d.Loads:
            name=load.Name().lower();bus=d.CktElement.BusNames()[0].split(".")[0].lower()
            load_nodes[name]=[f"{bus}.{n}" for n in d.CktElement.NodeOrder() if n];loads_at_bus[bus].append(name)
        static=dict(names=names,index=index,nt=nt,nc=nc,elem_start=elem_start,term_first=term_first,term_elem=term_elem,
            term_start=term_start,kind=kind,terminals=terminals,line_rating=line_rating,term_rating=term_rating,
            incidence=dict(incidence),load_nodes=load_nodes,loads_at_bus=dict(loads_at_bus),limits=None)
        self._static[region]=static
        return static

    def _flows(self, static):
        d=self.dss
        powers=np.asarray(d.PDElements.AllPowers(),dtype=float);currents=np.asarray(d.PDElements.AllCurrentsMagAng(),dtype=float)
        if len(powers)!=2*int(static["elem_start"][-1]+static["nt"][-1]*static["nc"][-1]) or len(currents)!=len(powers):
            raise ValueError("PD element batch arrays changed layout")
        if not np.all(np.isfinite(powers)) or not np.all(np.isfinite(currents)):raise PowerInfeasible("Nonfinite element flows")
        pterm=np.add.reduceat(powers[0::2],static["term_start"]);qterm=np.add.reduceat(powers[1::2],static["term_start"])
        return pterm,np.hypot(pterm,qterm),np.maximum.reduceat(currents[0::2],static["elem_start"])

    def _reference(self, region):
        settings={k:self.cfg[k] for k in ("voltage_min_pu","voltage_max_pu","voltage_tolerance","thermal_tolerance")}
        key=digest(dict(region=region,engine=self.version,settings=settings,limit_reference=self.cfg.get("limit_reference")))
        if key not in self._references:
            path=self.cache/region/f"reference-{key[:20]}.json"
            self._references[key]=cached_json(path,lambda:self._compute_reference(region,key))
        return self._references[key]

    def _compute_reference(self, region, key):
        """The published healthy base case at nominal load defines what is normal.

        Elements it already loads above nameplate, and load nodes it already
        holds outside the voltage band, keep their base-case value as limit.
        """
        if self.cfg.get("limit_reference")!="healthy_base_case":raise ValueError("Unknown power.limit_reference")
        d=self.dss;p=self.cfg
        with phase(f"AC {region}: healthy base-case reference at nominal load",enabled=self.verbose):
            self._compile(region,{},1.0,[]);d.Solution.Solve()
        if not d.Solution.Converged() or not d.Solution.ControlActionsDone():
            raise PowerInfeasible(f"{region}: healthy base case at nominal load does not converge")
        static=self._region_static(region);pterm,kva,current=self._flows(static)
        names=static["names"];kind=static["kind"];term_first=static["term_first"];nt=static["nt"]
        line_amps={};transformer_kva={};line_ratios=[];transformer_ratios=[]
        isolation={n.lower() for a in self.catalog["assets"].values() if a["kind"]=="power" and a["region"]==region for n in a["isolation_transformers"]}
        for i in np.flatnonzero(kind==1):
            rating=static["line_rating"][i];flow=float(current[i])
            if flow>(rating if rating>0 else 1e-6):line_amps[names[i]]=flow;line_ratios.append(flow/rating if rating>0 else math.inf)
        for i in np.flatnonzero(kind==2):
            span=range(int(term_first[i]),int(term_first[i]+nt[i]))
            over=[float(kva[j])>(static["term_rating"][j] if static["term_rating"][j]>0 else 1e-6) for j in span]
            if any(over) or names[i] in isolation:transformer_kva[names[i]]=[float(kva[j]) for j in span]
            if any(over):transformer_ratios.append(max(float(kva[j])/static["term_rating"][j] if static["term_rating"][j]>0 else math.inf for j in span))
        voltage=dict(zip((n.lower() for n in d.Circuit.AllNodeNames()),map(float,d.Circuit.AllBusMagPu())))
        low={};high={};dark=0
        for name,nodes in static["load_nodes"].items():
            volts=[voltage.get(n,0.0) for n in nodes]
            if not volts or not all(v>1e-4 for v in volts):dark+=1;continue
            for node,v in zip(nodes,volts):
                if v<p["voltage_min_pu"]:low[node]=v
                elif v>p["voltage_max_pu"]:high[node]=v
        finite=lambda xs:max((x for x in xs if math.isfinite(x)),default=None)
        return dict(digest=key,region=region,limit_reference="healthy_base_case",line_amps=line_amps,transformer_kva=transformer_kva,
            voltage_low=low,voltage_high=high,healthy_dark_loads=dark,
            preexisting=dict(lines_above_nameplate=len(line_amps),transformers_above_nameplate=len(transformer_ratios),
                max_line_loading_ratio=finite(line_ratios),max_transformer_loading_ratio=finite(transformer_ratios),
                load_nodes_below_voltage_min=len(low),load_nodes_above_voltage_max=len(high),
                min_load_voltage_pu=min(low.values(),default=None),max_load_voltage_pu=max(high.values(),default=None)))

    def _limits(self, region, reference):
        static=self._static[region]
        if static["limits"] is not None and static["limits"][0]==reference["digest"]:return static["limits"][1:]
        line=np.where(np.isfinite(static["line_rating"]),static["line_rating"],0.)
        term=np.where(np.isfinite(static["term_rating"]),static["term_rating"],0.)
        for name,amps in reference["line_amps"].items():
            i=static["index"][name];line[i]=max(line[i],amps)
        for name,values in reference["transformer_kva"].items():
            i=static["index"][name]
            for t,value in enumerate(values):j=static["term_first"][i]+t;term[j]=max(term[j],value)
        static["limits"]=(reference["digest"],line,term)
        return line,term

    def _measure_local(self, region, scale, closed, deratings, multipliers, reference):
        d=self.dss;p=self.cfg;static=self._static[region]
        voltage=dict(zip((n.lower() for n in d.Circuit.AllNodeNames()),map(float,d.Circuit.AllBusMagPu())))
        if not voltage or not all(math.isfinite(v) and v>=0 for v in voltage.values()):raise PowerInfeasible("Nonfinite voltages")
        tol=p["voltage_tolerance"];low=reference["voltage_low"];high=reference["voltage_high"]
        zones=defaultdict(float);subs=defaultdict(float);powered={};energized_loads=set();actual_kw={}
        actual_total=served=nominal=0.0;dark=0;energized_values=[];voltage_violations=[];missing=set(self.loads[region])
        for load in d.Loads:
            name=load.Name().lower();meta=self.loads[region].get(name)
            if meta is None:raise ValueError(f"Unexpected regional load {region}/{name}")
            missing.remove(name);kw=meta["kw"];m=multipliers.get(name,1.0);nominal+=kw
            if not math.isclose(load.kW(),kw*m,rel_tol=1e-8,abs_tol=1e-7) or not math.isclose(load.kvar(),float(meta["kvar"])*m,rel_tol=1e-8,abs_tol=1e-7):
                raise ValueError(f"Native P/Q differs from catalog times shedding multiplier: {region}/{name}")
            nodes=static["load_nodes"][name];volts=[voltage.get(n,0.0) for n in nodes]
            energized=bool(volts) and all(v>1e-4 for v in volts)
            power=d.CktElement.Powers();nc=d.CktElement.NumConductors()
            actual=sum(float(x) for x in power[:2*nc:2])
            if not math.isfinite(actual):raise PowerInfeasible("Nonfinite load power")
            actual_total+=actual
            delivered=max(0.0,min(kw,actual)) if energized else 0.0
            served+=delivered;zones[meta["zone"]]+=delivered;subs[meta["substation_id"]]+=delivered
            dark+=not energized;in_band=energized
            if energized:
                energized_loads.add(name);actual_kw[name]=actual;energized_values.extend(volts)
                for node,v in zip(nodes,volts):
                    lower=min(p["voltage_min_pu"],low.get(node,math.inf))-tol;upper=max(p["voltage_max_pu"],high.get(node,-math.inf))+tol
                    if v<lower:voltage_violations.append(("low",name,node,v,lower));in_band=False
                    elif v>upper:voltage_violations.append(("high",name,node,v,upper));in_band=False
            good=in_band and delivered/max(kw,1e-12)>=p["signal_min_served_fraction"]
            for sid in self.signal_loads.get((region,name),[]):powered[sid]=bool(good)
        if missing:raise ValueError("Native load identities missing from solve")
        pterm,kva,current=self._flows(static);line_normal,term_normal=self._limits(region,reference)
        emergency=p["emergency_factor"];line_limit=line_normal*emergency;term_limit=term_normal*emergency
        for name,factor in deratings.items():
            # A damaged unit keeps only its remaining share of the normal rating; no emergency loading.
            i=static["index"][name.lower()];a=static["term_first"][i];b=a+static["nt"][i];term_limit[a:b]=term_normal[a:b]*factor
        kind=static["kind"];names=static["names"];thermal=[];unrated=[];max_loading=0.0
        ratio=np.zeros(len(names))
        lines=kind==1
        with np.errstate(divide="ignore",invalid="ignore"):
            line_ratio=np.where(line_limit>0,current/line_limit,np.where(current>1e-6,np.inf,0.))
            term_ratio=np.where(term_limit>0,kva/term_limit,np.where(kva>1e-6,np.inf,0.))
        ratio[lines]=line_ratio[lines]
        transformer_terms=kind[static["term_elem"]]==2
        np.maximum.at(ratio,static["term_elem"][transformer_terms],term_ratio[transformer_terms])
        finite=ratio[np.isfinite(ratio)]
        max_loading=float(finite.max()) if finite.size else 0.0
        for i in np.flatnonzero(ratio>1+p["thermal_tolerance"]):
            a=int(static["term_first"][i]);span=range(a,a+int(static["nt"][i]))
            entering=max(float(pterm[j]) for j in span)
            exits=[static["terminals"][i][j-a] for j in span if pterm[j]< -1e-3]
            if not math.isfinite(ratio[i]):unrated.append(names[i])
            thermal.append((int(i),float(ratio[i]),entering,exits))
        source_kw=-float(d.Circuit.TotalPower()[0]);losses_kw=float(d.Circuit.Losses()[0])/1000
        residual=source_kw-actual_total-losses_kw
        balance_ok=abs(residual)<=max(.1,abs(source_kw)*p["balance_relative_tolerance"])
        feasible=balance_ok and not thermal and not voltage_violations
        result=dict(region=region,feasible=bool(feasible),converged=True,controls_settled=True,
            load_scale=scale,served_kw=served,nominal_kw=nominal,actual_load_kw=actual_total,source_kw=source_kw,
            losses_kw=losses_kw,active_balance_residual_kw=residual,
            min_loaded_voltage_pu=min(energized_values) if energized_values else None,
            max_loaded_voltage_pu=max(energized_values) if energized_values else None,
            max_thermal_loading_ratio=max_loading,unrated_energized_elements=unrated,dark_loads=dark,
            thermal_violations=len(thermal),voltage_violations=len(voltage_violations),
            zone_served_kw=dict(zones),substation_served_kw=dict(subs),signal_powered=powered,closed_switches=list(closed),
            voltage_limits=[p["voltage_min_pu"],p["voltage_max_pu"]],engine=self.version)
        return dict(result=result,thermal=thermal,voltage=voltage_violations,energized=energized_loads,actual_kw=actual_kw,pterm=pterm)

    def _local_solve(self, region, damage, scale, closed, multipliers, reference):
        deratings,disabled=self._compile(region,damage,scale,closed,multipliers)
        d=self.dss;started=time.monotonic()
        self._mark("solve",region=region,scale=scale,damage=damage,closed=closed,shed_loads=len(multipliers))
        with phase(f"AC {region}: power flow; local shedding; ties={len(closed)}",enabled=self.verbose):d.Solution.Solve()
        self.stats["solve_seconds"]+=time.monotonic()-started;self.stats["solve_calls"]+=1
        if not d.Solution.Converged() or not d.Solution.ControlActionsDone():
            raise PowerInfeasible("Power flow or controller iteration did not converge")
        for name in disabled:
            d.Circuit.SetActiveElement("Transformer."+name)
            if d.CktElement.Name().lower()!="transformer."+name or d.CktElement.Enabled():raise ValueError("Fault isolation failed: "+name)
        started=time.monotonic();self._mark("measure",region=region,scale=scale,damage=damage,closed=closed)
        measured=self._measure_local(region,scale,closed,deratings,multipliers,reference)
        self.stats["measure_seconds"]+=time.monotonic()-started
        measured["result"].update(isolated_transformers=sorted(disabled),rating_derating_factors=deratings)
        return measured

    def _shed(self, region, measured, multipliers, pickup):
        """Lower multipliers below violated elements; False when shedding cannot help.

        With ``pickup`` (a tie trial) only loads newly energized by the tie may
        be reduced; customers supplied before the tie are never shed for it.
        """
        p=self.cfg;static=self._static[region];target={}
        if any(v[0]=="high" for v in measured["voltage"]):return False
        def lower(name,value):
            if value<min(target.get(name,math.inf),multipliers.get(name,1.0))-1e-12:target[name]=max(0.0,value)
        for element,ratio,entering,exits in measured["thermal"]:
            loads=[n for n in downstream_loads(exits,static["incidence"],measured["pterm"],static["term_first"],static["terminals"],
                static["loads_at_bus"],element) if n in measured["energized"]]
            if pickup is not None:
                loads=[n for n in loads if n in pickup]
                if not loads:return False
                needed=entering*(1-p["shed_safety"]/ratio) if math.isfinite(ratio) else entering
                available=sum(max(0.0,measured["actual_kw"][n]) for n in loads)
                factor=0.0 if available<=needed else 1-needed/available
            else:
                if not loads:return False
                factor=p["shed_safety"]/ratio if math.isfinite(ratio) else 0.0
            for name in loads:lower(name,multipliers.get(name,1.0)*factor)
        lows=[v for v in measured["voltage"] if v[0]=="low"]
        if lows:
            if pickup is not None:
                if not pickup:return False
                for name in pickup:lower(name,multipliers.get(name,1.0)*p["voltage_shed_step"])
            else:
                feeders={self.loads[region][v[1]]["feeder_id"] for v in lows}
                for name,meta in self.loads[region].items():
                    if meta["feeder_id"] in feeders and name in measured["energized"]:lower(name,multipliers.get(name,1.0)*p["voltage_shed_step"])
        multipliers.update(target)
        return bool(target)

    def _local_state(self, region, damage, closed, reference, incumbent=None, scale=1.0):
        p=self.cfg;attempts=[]
        multipliers=dict(incumbent["_multipliers"]) if incumbent is not None else {}
        protected=incumbent["_energized"] if incumbent is not None else None;pickup=None
        for iteration in range(p["max_shed_iterations"]+1):
            try:measured=self._local_solve(region,damage,scale,closed,multipliers,reference)
            except (PowerInfeasible,self.dss.DSSException) as exc:
                attempts.append(dict(iteration=iteration,error=str(exc)));break
            result=measured["result"]
            if protected is not None and pickup is None:pickup=measured["energized"]-protected
            attempts.append(dict(iteration=iteration,load_scale=scale,feasible=result["feasible"],thermal_violations=result["thermal_violations"],
                voltage_violations=result["voltage_violations"],max_thermal_loading_ratio=result["max_thermal_loading_ratio"],
                shed_loads=sum(1 for v in multipliers.values() if v<1)))
            if result["feasible"]:
                result.update(curtailment_trials=attempts,shed_iterations=iteration,regional_fallback=False,
                              _multipliers=dict(multipliers),_energized=measured["energized"])
                return result
            if iteration==p["max_shed_iterations"] or not self._shed(region,measured,multipliers,pickup):break
        if incumbent is not None:
            raise PowerInfeasible(f"{region}: tie pickup does not fit without shedding supplied customers; trials={attempts}")
        # Last resort, recorded: the legacy uniform grid on top of the local multipliers.
        for grid in self.cfg["load_scales"]:
            if grid>=scale:continue
            try:measured=self._local_solve(region,damage,grid,closed,multipliers,reference)
            except (PowerInfeasible,self.dss.DSSException) as exc:
                attempts.append(dict(load_scale=grid,error=str(exc)));continue
            result=measured["result"]
            attempts.append(dict(load_scale=grid,feasible=result["feasible"],thermal_violations=result["thermal_violations"],
                voltage_violations=result["voltage_violations"],max_thermal_loading_ratio=result["max_thermal_loading_ratio"]))
            if result["feasible"]:
                self.stats["regional_fallbacks"]+=1
                result.update(curtailment_trials=attempts,shed_iterations=len(attempts),regional_fallback=True,
                              _multipliers=dict(multipliers),_energized=measured["energized"])
                return result
        raise PowerInfeasible(f"{region}: no feasible operating point with local shedding or regional fallback; trials={attempts}")

    def _operate_local(self, region, damage):
        started=time.monotonic();before=dict(self.stats)
        reference=self._reference(region)
        best=self._local_state(region,damage,[],reference);evaluations=[]
        for iteration in range(self.cfg["max_switch_closures"]):
            # Reestablish the best state before identifying energized/dark endpoints.
            if iteration:self._local_solve(region,damage,best["load_scale"],best["closed_switches"],best["_multipliers"],reference)
            candidates=self._tie_candidates(best["closed_switches"])[:self.cfg["max_switch_trials"]]
            incumbent=best
            for name in candidates:
                try:
                    trial=self._local_state(region,damage,best["closed_switches"]+[name],reference,incumbent=best,scale=best["load_scale"])
                    evaluations.append(dict(switch=name,served_kw=trial["served_kw"],feasible=True))
                    if trial["served_kw"]>incumbent["served_kw"]+1e-6:incumbent=trial
                except PowerInfeasible as exc:evaluations.append(dict(switch=name,feasible=False,error=str(exc)))
            if incumbent is best:break
            best=incumbent
        multipliers=best.pop("_multipliers");best.pop("_energized")
        meta=self.loads[region]
        best.update(switch_trials=evaluations,damage=damage,limit_reference=reference["digest"],
            preexisting_base_case=reference["preexisting"],shed_loads=sum(1 for v in multipliers.values() if v<1-1e-12),
            shed_nominal_kw=sum(meta[n]["kw"]*(1-v) for n,v in multipliers.items()),
            operator="healthy-base-case-referenced limits; downstream load shedding; tie pickup only within spare capacity; no optimality claim")
        best["performance"]={k:v-before.get(k,0) for k,v in self.stats.items()}
        best["performance"]["wall_seconds"]=time.monotonic()-started
        self._mark("complete",region=region,damage=damage,performance=best["performance"])
        return best
