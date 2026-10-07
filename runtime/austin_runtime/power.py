"""Physical three-phase OpenDSS operation on TAMU's complete regional circuits.

The operator searches an explicit uniform regional curtailment grid and existing
open switches. It is a feasible heuristic, not an OPF or guaranteed optimum.
"""
from __future__ import annotations
import concurrent.futures as futures
import math, multiprocessing, os, signal, subprocess, sys, tempfile, time
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
        self._changed_enabled={}; self._changed_terminals={}
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

    def _compile(self, region, damage, scale, closed):
        self._load_circuit(region,scale,closed)
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
