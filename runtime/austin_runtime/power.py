"""Physical three-phase OpenDSS operation on TAMU's complete regional circuits.

The operator searches an explicit uniform regional curtailment grid and existing
open switches. It is a feasible heuristic, not an OPF or guaranteed optimum.
"""
from __future__ import annotations
import math
from collections import defaultdict
from pathlib import Path
from .common import cached_json, digest, local_path, rows


class PowerInfeasible(RuntimeError):
    pass


class PowerEngine:
    def __init__(self, cfg, catalog, fingerprint, scratch):
        from opendssdirect import dss
        self.cfg=cfg["power"]; self.catalog=catalog; self.fingerprint=fingerprint
        self.prepared=local_path(cfg,"prepared"); self.cache=local_path(cfg,"cache")/fingerprint/"power"
        self.scratch=Path(scratch); self.scratch.mkdir(parents=True,exist_ok=True)
        self.dss=dss.NewContext()
        self.version=self.dss.Basic.Version()
        self.loads=defaultdict(dict)
        for r in rows(self.prepared/"load_catalog.csv.gz"):
            r["kw"]=float(r["kw"]); self.loads[r["region"]][r["load_id"]]=r
        self.signal_loads=defaultdict(list)
        for r in catalog["signals"]:
            region=catalog["substation_region"][r["substation_id"]]
            self.signal_loads[(region,r["load_id"].lower())].append(r["signal_id"])

    def evaluate(self, remaining):
        results=[]
        for region in self.catalog["regions"]:
            damage={key:factor for key,factor in remaining.items() if key.startswith("power:") and self.catalog["assets"][key]["region"]==region}
            key=digest(dict(region=region,damage=damage,engine=self.version,settings=self.cfg))
            results.append(cached_json(self.cache/region/(key+".json"),lambda r=region,d=damage:self._operate(r,d)))
        zones=defaultdict(float); subs=defaultdict(float); signals={}
        for r in results:
            for k,v in r["zone_served_kw"].items():zones[k]+=v
            for k,v in r["substation_served_kw"].items():subs[k]+=v
            signals.update(r["signal_powered"])
        return dict(served_kw=sum(r["served_kw"] for r in results),nominal_kw=sum(r["nominal_kw"] for r in results),
            zone_served_kw=dict(zones),substation_served_kw=dict(subs),signal_powered=signals,
            regions=[{k:v for k,v in r.items() if k not in ("zone_served_kw","substation_served_kw","signal_powered")} for r in results])

    def _compile(self, region, damage, scale, closed):
        d=self.dss
        d.Basic.AllowChangeDir(False)
        d.Basic.DataPath(str(self.scratch))
        master=self.prepared/"models"/region/"Master_runtime.dss"
        d(f'Redirect "{master}"')
        d.Basic.DataPath(str(self.scratch))
        d("Set mode=snapshot controlmode=static maxcontroliter=100 maxiterations=100")
        d.Solution.LoadMult(scale)
        disabled=set(); deratings={}
        for key,factor in damage.items():
            for name in self.catalog["assets"][key]["isolation_transformers"]:
                if factor==0:
                    d("Disable "+name); disabled.add(name.split(".",1)[1])
                else: deratings[name]=factor
        # Only disable a regulator controlling an isolated transformer. Cross-feeder
        # references elsewhere remain intact; no control object is deleted.
        for name in d.RegControls.AllNames():
            if name.lower()=="none":continue
            d.RegControls.Name(name)
            if d.RegControls.Transformer().lower().removeprefix("transformer.") in disabled: d("Disable RegControl."+name)
        for name in closed:
            d.Circuit.SetActiveElement(name)
            if d.CktElement.Name().lower()!=name.lower():raise ValueError("Switch not found: "+name)
            if not d.CktElement.Enabled():raise ValueError("Cannot close disabled component: "+name)
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
                self.dss.Solution.Solve()
                if not self.dss.Solution.Converged() or not self.dss.Solution.ControlActionsDone():
                    raise PowerInfeasible("Power flow or controller iteration did not converge")
                for name in disabled:
                    self.dss.Circuit.SetActiveElement("Transformer."+name)
                    if self.dss.CktElement.Name().lower()!="transformer."+name or self.dss.CktElement.Enabled():
                        raise ValueError("Fault isolation failed: "+name)
                result=self._measure(region,scale,closed,deratings)
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
        return best
