from __future__ import annotations
import gzip, heapq, json, math, os, re, shutil, subprocess, tempfile, time
from collections import OrderedDict, defaultdict
from pathlib import Path
import numpy as np
from .common import AUSTIN, RUNTIME, atomic_json, digest, file_lock, local_path, rows, sha_file
from .progress import phase


class TrafficState:
    def __init__(self, links, costs, flows, factors, report, first_thru=1118):
        self.costs=costs; self.flows=flows; self.factors=factors; self.report=report; self.first_thru=first_thru
        self.forward=defaultdict(list); self.reverse=defaultdict(list); self._distances={}
        for r,c,f in zip(links,costs,factors):
            if f<=0:continue # Closed arcs remain a penalty device in assignment, never in crew routing.
            u,v,i=int(r["from_node"]),int(r["to_node"]),int(r["link_id"])
            self.forward[u].append((v,float(c),i)); self.reverse[v].append((u,float(c),i))

    def distances(self, sources, reverse=False):
        sources=tuple(sorted(set(sources))); key=(sources,reverse)
        if key in self._distances:return self._distances[key]
        adjacency=self.reverse if reverse else self.forward
        dist={s:0.0 for s in sources}; heap=[(0.0,s) for s in sources]; heapq.heapify(heap)
        while heap:
            value,u=heapq.heappop(heap)
            if value!=dist[u]:continue
            # A centroid may be the origin/destination; it can never be an intermediate road node.
            if u<self.first_thru and u not in sources:continue
            for v,c,_ in adjacency.get(u,()):
                new=value+c
                if new<dist.get(v,math.inf):dist[v]=new;heapq.heappush(heap,(new,v))
        self._distances[key]=dist
        return dist

    def travel(self, origin, targets):
        dist=self.distances([origin]); reachable=[(dist.get(t,math.inf),t) for t in targets]
        return min(reachable) if reachable else (math.inf,origin)


class TrafficEngine:
    def __init__(self,cfg,catalog,fingerprint,scratch,*,verbose=False):
        self.verbose=verbose
        self.cfg=cfg; self.catalog=catalog; self.scratch=Path(scratch);self.scratch.mkdir(parents=True,exist_ok=True)
        self.links=list(rows(AUSTIN/"data/processed/road/links.csv"))
        self.net=AUSTIN/"data/processed/road/Austin_net.tntp";self.trips=AUSTIN/"data/processed/road/Austin_trips.tntp"
        self.binary=RUNTIME/"build/tap-b/bin/tap"
        if not self.binary.is_file():raise FileNotFoundError("Build native TAP-B with runtime/setup.sh on the destination host")
        self.binary_hash=sha_file(self.binary)
        self.cache=local_path(cfg,"cache")/fingerprint/"traffic"; self.memory=OrderedDict()
        self.supply=np.zeros(7467,dtype=float)
        for row in rows(AUSTIN/"data/processed/road/od.csv.gz"):
            o,d,q=int(row["origin"]),int(row["destination"]),float(row["demand_source"])
            self.supply[o]+=q;self.supply[d]-=q
        self.signal_links=defaultdict(set)
        for r in catalog["approaches"]:self.signal_links[r["signal_id"]].add(int(r["link_id"])-1)

    def evaluate(self,remaining,signal_powered,penalty=None):
        factors=np.ones(len(self.links),dtype=float)
        for key,factor in remaining.items():
            if key.startswith("road:"):
                for link in self.catalog["assets"][key]["link_ids"]:factors[link-1]*=factor
        off_links=set()
        for sid,indices in self.signal_links.items():
            if not signal_powered.get(sid,False):off_links.update(indices)
        for i in off_links:factors[i]*=self.cfg["traffic"]["signal_outage_capacity_factor"]
        penalty=penalty or self.cfg["traffic"]["closure_penalty"]
        key=digest(dict(factors=factors.tolist(),penalty=penalty,binary=self.binary_hash,gap=self.cfg["traffic"]["relative_gap"]))
        if key in self.memory:self.memory.move_to_end(key);return self.memory[key]
        path=self.cache/(key+".npz");path.parent.mkdir(parents=True,exist_ok=True)
        with file_lock(path.with_suffix(".lock")):
            if path.exists():
                with np.load(path,allow_pickle=False) as data:
                    costs=data["costs"].copy();flows=data["flows"].copy();report=json.loads(str(data["report"]))
            else:
                costs,flows,report=self._solve(factors,penalty,key)
                fd,name=tempfile.mkstemp(dir=path.parent,prefix=key+".",suffix=".tmp")
                try:
                    with os.fdopen(fd,"wb") as f:
                        np.savez_compressed(f,costs=costs,flows=flows,report=json.dumps(report,allow_nan=False))
                        f.flush();os.fsync(f.fileno())
                    os.replace(name,path)
                finally:
                    if os.path.exists(name):os.unlink(name)
        state=TrafficState(self.links,costs,flows,factors,report)
        self.memory[key]=state
        if len(self.memory)>8:self.memory.popitem(last=False)
        return state

    def _solve(self,factors,penalty,key):
        # Upstream TAP-B writes flows.txt. Retain ALL files if solving or parsing fails.
        folder=Path(tempfile.mkdtemp(prefix="tapb-",dir=self.scratch))
        try:
            network=folder/"network.tntp"
            head=self.net.read_text().splitlines()[:8]
            lines=[]
            for r,factor in zip(self.links,factors):
                values=[r[k] for k in ("from_node","to_node","capacity_source","length_source","free_flow_time_source","bpr_alpha","bpr_beta","speed_source","toll_source","link_type")]
                if factor<=0:
                    values[2]="1000000000";values[4]=str(penalty);values[5]="0"
                else:values[2]=format(float(values[2])*factor,".15g")
                lines.append("\t".join(values)+"\t;")
            network.write_text("\n".join(head+lines)+"\n")
            args=[str(self.binary),str(self.cfg["traffic"]["relative_gap"]),"1",str(network),str(self.trips),str(self.cfg["runtime"]["tapb_threads"])]
            start=time.monotonic();logfile=folder/"tapb.log"
            try:
                with phase(f"TAP-B solve; log={logfile}",enabled=self.verbose), logfile.open("w") as handle:
                    result=subprocess.run(args,cwd=folder,stdout=handle,stderr=subprocess.STDOUT,
                        timeout=self.cfg["traffic"]["timeout_seconds"],check=False)
            except subprocess.TimeoutExpired:
                error=self.scratch/(key+".timeout.log");error.write_text(logfile.read_text())
                raise RuntimeError(f"TAP-B timed out; {error}") from None
            log=logfile.read_text();gaps=re.findall(r"Iteration\s+(\d+):\s+gap\s+([\d.eE+-]+)",log)
            gap=float(gaps[-1][1]) if gaps else math.inf
            if result.returncode or not math.isfinite(gap) or not 0<=gap<=self.cfg["traffic"]["relative_gap"]:
                error=self.scratch/(key+".failed.log");error.write_text(log)
                raise RuntimeError(f"TAP-B failed convergence/exit check; {error}")
            flow_file=folder/"flows.txt"
            if not flow_file.is_file():
                raise RuntimeError(f"TAP-B returned success but did not write flows.txt; files={sorted(p.name for p in folder.iterdir())}")
            records=[]
            for line in flow_file.read_text().splitlines():
                m=re.match(r"\((\d+),(\d+)\)\s+(\S+)\s+(\S+)",line.strip())
                if m:records.append(tuple(map(float,m.groups())))
            if len(records)!=len(self.links):raise ValueError("TAP-B lost arcs")
            costs=[];flows=[];balance=np.zeros_like(self.supply)
            for row,(u,v,flow,cost) in zip(self.links,records):
                if (int(u),int(v))!=(int(row["from_node"]),int(row["to_node"])):raise ValueError("TAP-B changed arc order")
                if not math.isfinite(flow+cost) or flow<-.000001 or cost<0:raise ValueError("Invalid TAP-B flow/cost")
                costs.append(cost);flows.append(flow);balance[int(u)]+=flow;balance[int(v)]-=flow
            error=float(np.max(np.abs(balance-self.supply)))
            if error>1e-4:raise ValueError(f"TAP-B flow conservation failed: {error}")
            tstt=float(re.findall(r"Aggregate TSTT:\s+([\d.eE+-]+)",log)[-1])
            if not math.isclose(float(np.dot(costs,flows)),tstt,rel_tol=1e-6):raise ValueError("TSTT mismatch")
            report=dict(tstt=tstt,relative_gap=gap,iterations=int(gaps[-1][0]),wall_seconds=time.monotonic()-start,
                max_node_imbalance=error,parallel_arcs_preserved=True,network_sha256=sha_file(network),binary_sha256=self.binary_hash,
                closure_penalty=penalty,closed_arc_flow_sum=float(np.sum(np.asarray(flows)[factors<=0])),
                closure_semantics="finite assignment penalty; closed arcs excluded from crew routes; flow sum is NOT unique unserved OD demand")
            answer=np.asarray(costs),np.asarray(flows),report
        except BaseException as exc:
            (folder/"failure.txt").write_text(f"{type(exc).__name__}: {exc}\n")
            # Do not let TemporaryDirectory erase the evidence of an interface failure.
            if isinstance(exc,Exception):
                raise RuntimeError(f"{exc}; TAP-B diagnostic files retained at {folder}") from exc
            raise
        else:
            shutil.rmtree(folder)
            return answer
