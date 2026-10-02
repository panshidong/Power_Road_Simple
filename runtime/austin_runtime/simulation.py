from __future__ import annotations
import math
from pathlib import Path
import numpy as np
from .power import PowerEngine
from .traffic import TrafficEngine


def ratio(reference, current):
    if not math.isfinite(reference) or not math.isfinite(current):return 0.0
    if reference==current==0:return 1.0
    return max(0.0,min(1.0,reference/max(current,1e-12)))


def gini(values):
    x=np.sort(np.asarray(values,dtype=float)); n=len(x)
    return float(2*np.dot(np.arange(1,n+1),x)/(n*x.sum())-(n+1)/n) if n and x.sum()>0 else 0.0


class Coupled:
    def __init__(self,cfg,catalog,fingerprint,scratch,*,power_workers=1,verbose=False):
        self.cfg=cfg;self.catalog=catalog
        self.power=PowerEngine(cfg,catalog,fingerprint,Path(scratch)/"dss",region_workers=power_workers,verbose=verbose)
        try:
            self.traffic=TrafficEngine(cfg,catalog,fingerprint,Path(scratch)/"tapb",verbose=verbose)
            self.healthy_power=self.power.evaluate({})
            if self.healthy_power["served_kw"]<=0:raise RuntimeError("Healthy AC case has no served load")
            self.healthy_traffic=self.traffic.evaluate({},self.healthy_power["signal_powered"])
            self.depot=catalog["depot"];self.shelter=catalog["shelter"]
            self.essential=[self.shelter]+[catalog["assets"]["power:"+s]["targets"][0] for s in catalog["critical_substations"]]
            self.baseline_tt=self.healthy_traffic.distances(self.essential,reverse=True)
            self.zones=sorted(int(z) for z,v in self.healthy_power["zone_served_kw"].items() if v>1e-9 and math.isfinite(self.baseline_tt.get(int(z),math.inf)))
            if not self.zones:raise RuntimeError("No zones have both baseline power and essential access")
            self.excluded_zones=sorted(set(map(int,catalog["zone_nominal_kw"]))-set(self.zones))
        except BaseException:
            self.power.close(cancel=True)
            raise

    def close(self, *, cancel=False):
        self.power.close(cancel=cancel)

    def state(self,remaining,penalty=None):
        power=self.power.evaluate(remaining)
        traffic=self.traffic.evaluate(remaining,power["signal_powered"],penalty)
        p=max(0.0,min(1.0,power["served_kw"]/self.healthy_power["served_kw"]))
        r=ratio(self.healthy_traffic.report["tstt"],traffic.report["tstt"])
        tt=traffic.distances(self.essential,reverse=True)
        electric=[max(0.,min(1.,power["zone_served_kw"].get(str(z),0)/self.healthy_power["zone_served_kw"][str(z)])) for z in self.zones]
        access=[ratio(self.baseline_tt[z],tt.get(z,math.inf)) for z in self.zones]
        return dict(power=power,traffic=traffic,power_func=p,road_func=r,electric=electric,access=access)

    def site_travel(self,traffic,origin,asset,roundtrip_lastmile=False):
        info=self.catalog["assets"][asset];cost,target=traffic.travel(origin,info["targets"])
        extra=info["offroad_m"]/(self.cfg["recovery"]["offroad_speed_kph"]*1000/60)
        return cost*self.cfg["traffic"]["time_to_minutes"]+extra*(2 if roundtrip_lastmile else 1),target

    def coalition(self,scenario,repaired):
        remaining={k:v for k,v in scenario["damage"].items() if k not in repaired}
        state=self.state(remaining)
        access=[]
        for asset in scenario["damage"]:
            if not asset.startswith("power:"):continue
            base,_=self.site_travel(self.healthy_traffic,self.depot,asset)
            current,_=self.site_travel(state["traffic"],self.depot,asset)
            access.append(ratio(base,current))
        return .5*(state["power_func"]+state["road_func"]),sum(access)/len(access) if access else 1.0

    def record(self,time,remaining,state,critical):
        power=state["power"];multiplier=self.cfg["recovery"]["critical_load_multiplier"]
        numerator=denominator=0.0
        for sid,base in self.healthy_power["substation_served_kw"].items():
            weight=multiplier if sid in critical else 1.0
            denominator+=weight*base;numerator+=weight*min(base,power["substation_served_kw"].get(sid,0))
        reference,_=self.healthy_traffic.travel(self.depot,[self.shelter]);current,_=state["traffic"].travel(self.depot,[self.shelter])
        return dict(time=time,remaining=dict(remaining),power_func=state["power_func"],road_func=state["road_func"],
            power_absolute=power["served_kw"]/power["nominal_kw"],served_kw=power["served_kw"],
            weighted_power_func=numerator/max(denominator,1e-12),shelter_access=ratio(reference,current),
            electric=state["electric"],access=state["access"],traffic_quality=state["traffic"].report,
            power_operation=power["regions"])

    def simulate(self,scenario,sequence,*,greedy=False,penalty=None,critical=None,equity=None):
        if set(sequence)!=set(scenario["damage"]) or len(sequence)!=len(set(sequence)):
            raise ValueError("Recovery sequence must contain each initially damaged asset once")
        recovery=self.cfg["recovery"];remaining=dict(scenario["damage"]);pending=list(sequence)
        critical=critical if critical is not None else self.catalog["critical_substations"]
        crews=[]
        for kind in ("power","road"):
            for i in range(recovery[kind+"_crews"]):
                crews.append(dict(id=f"{kind}-{i}",kind=kind,node=self.depot,job=None,finish=None,target=None))
        time=0.0;events=[];dispatch=[];state=self.state(remaining,penalty);stop="completed"
        while True:
            events.append(self.record(time,remaining,state,critical))
            if not remaining:break
            # Repairs in progress stay failed until their completion event.
            for crew in crews:
                if crew["job"] is not None:continue
                eligible=[]
                for rank,asset in enumerate(pending):
                    if self.catalog["assets"][asset]["kind"]!=crew["kind"]:continue
                    travel,target=self.site_travel(state["traffic"],crew["node"],asset,roundtrip_lastmile=True)
                    if math.isfinite(travel):eligible.append((rank,asset,travel,target))
                if not eligible:continue
                gains={}
                if greedy:
                    current=.5*(state["power_func"]+state["road_func"])
                    for _,asset,_,_ in eligible:
                        trial=self.state({k:v for k,v in remaining.items() if k!=asset},penalty)
                        gains[asset]=.5*(trial["power_func"]+trial["road_func"])-current
                    choice=min(eligible,key=lambda item:(-gains[item[1]],item[1]))
                else:choice=eligible[0]
                _,asset,travel,target=choice;pending.remove(asset)
                finish=time+travel+recovery[crew["kind"]+"_repair_minutes"]
                dispatch.append(dict(time=time,crew=crew["id"],asset=asset,origin=crew["node"],target=target,
                    travel_minutes=travel,finish=finish,gain=gains.get(asset),candidate_gains=gains))
                crew.update(job=asset,finish=finish,target=target)
            busy=[c for c in crews if c["job"] is not None]
            if not busy:
                stop="unreachable_pending_jobs";time=recovery["horizon_minutes"]
                events.append(self.record(time,remaining,state,critical));break
            next_time=min(c["finish"] for c in busy)
            if next_time>recovery["horizon_minutes"]:
                stop="horizon";time=recovery["horizon_minutes"]
                events.append(self.record(time,remaining,state,critical));break
            time=next_time
            for crew in busy:
                if abs(crew["finish"]-time)<=1e-9:
                    del remaining[crew["job"]]
                    crew.update(node=crew["target"],job=None,finish=None,target=None)
            state=self.state(remaining,penalty)
        metrics=measure(events,self.zones,equity or self.cfg["equity"])
        return dict(scenario_id=scenario["id"],scenario_seed=scenario["seed"],sequence=sequence,events=events,
            dispatch=dispatch,metrics=metrics,stop_reason=stop,complete=not remaining,remaining=list(remaining),
            zones=self.zones,excluded_baseline_zones=self.excluded_zones,
            healthy_served_fraction=self.healthy_power["served_kw"]/self.healthy_power["nominal_kw"])


def measure(events,zones,eq):
    times=np.asarray([e["time"] for e in events]);dt=np.diff(times);end=float(times[-1]);n=len(zones)
    p=np.asarray([e["power_func"] for e in events]);r=np.asarray([e["road_func"] for e in events])
    wp=np.asarray([e["weighted_power_func"] for e in events]);shelter=np.asarray([e["shelter_access"] for e in events])
    electric=np.asarray([e["electric"] for e in events]);access=np.asarray([e["access"] for e in events])
    cri=eq["electric_weight"]*electric+eq["access_weight"]*access
    def first_hit(values,threshold):
        result=np.full(n,end);censored=[]
        for j,z in enumerate(zones):
            hits=np.flatnonzero(values[:,j]>=threshold)
            if len(hits):result[j]=times[hits[0]]
            else:censored.append(z)
        return result,censored
    restore,censored=first_hit(cri,eq["cri_threshold"])
    ar,acensored=first_hit(access,eq["access_threshold"])
    avg=np.sum(cri[:-1]*dt[:,None],axis=0)/end if end else cri[-1]
    result=dict(triangle_area=float(np.dot((1-p[:-1])+(1-r[:-1]),dt)),
        power_loss=float(np.dot(1-p[:-1],dt)),road_loss=float(np.dot(1-r[:-1],dt)),
        weighted_triangle_area=float(np.dot((1-r[:-1])+(1-wp[:-1])+(1-shelter[:-1]),dt)),
        gini_restore=gini(restore),var_restore=float(np.var(restore)),p90_restore=float(np.quantile(restore,.9)),
        p95_restore=float(np.quantile(restore,.95)),min_time_avg_cri=float(min(avg)),
        maximin_time_avg_cri_loss=float(1-min(avg)),mean_time_avg_cri=float(np.mean(avg)),
        p90_access_restore=float(np.quantile(ar,.9)),p95_access_restore=float(np.quantile(ar,.95)),
        gini_access_restore=gini(ar),gini_final_cri=gini(cri[-1]),mean_final_cri=float(np.mean(cri[-1])),
        min_final_cri=float(min(cri[-1])),shelter_access_initial=float(shelter[0]),shelter_access_final=float(shelter[-1]),
        completion_or_censor_time=end,censored_cri_zones=len(censored),censored_access_zones=len(acensored))
    if not all(math.isfinite(x) for x in result.values()):raise ValueError("Nonfinite recovery metric")
    return result
