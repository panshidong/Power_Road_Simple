from __future__ import annotations
import itertools, math, random
from collections import defaultdict
import networkx as nx
from .common import AUSTIN, atomic_json, cached_json, digest, local_path, read_json, rows


def construct(engine,scenario,cfg,checkpoint_path,signature):
    assets=sorted(scenario["damage"]); count=cfg["task_a"]["shapley_permutations"]
    progress=dict(signature=signature,permutations=0,sums={a:[0.,0.] for a in assets},squares={a:[0.,0.] for a in assets})
    if checkpoint_path.exists():
        progress=read_json(checkpoint_path)
        if progress["signature"]!=signature:raise ValueError("Construction checkpoint fingerprint mismatch")
    for index in range(progress["permutations"],count):
        permutation=assets.copy();random.Random(scenario["seed"]+303+index*1000003).shuffle(permutation)
        repaired=set(); previous=engine.coalition(scenario,repaired);initial=previous; total=[0.,0.]
        for a in permutation:
            repaired.add(a);current=engine.coalition(scenario,repaired)
            for j in (0,1):
                delta=current[j]-previous[j];progress["sums"][a][j]+=delta;progress["squares"][a][j]+=delta*delta;total[j]+=delta
            previous=current
        if any(abs(total[j]-(previous[j]-initial[j]))>1e-8 for j in (0,1)):raise ValueError("Shapley efficiency check failed")
        progress["permutations"]=index+1;atomic_json(checkpoint_path,progress)
    averages={a:[x/count for x in progress["sums"][a]] for a in assets}
    scores={"CEN":{a:engine.catalog["centrality"][a] for a in assets},"JSH":{a:v[0] for a,v in averages.items()}}
    for alpha in sorted(set(cfg["task_a"]["alpha_values"]+[1.0])):
        key="IJSH" if alpha==1 else f"IJSH_a{alpha:g}"
        scores[key]={a:(v[0]+alpha*v[1])/(1+alpha) for a,v in averages.items()}
    se={a:[math.sqrt(max(0.,(progress["squares"][a][j]-progress["sums"][a][j]**2/count)/max(count-1,1))/count) for j in (0,1)] for a in assets}
    return dict(scenario=scenario,scores=scores,permutations=count,mc_standard_errors_joint_and_access=se,
        method="E[scenario permutation marginal | asset damaged]; JSH and IJSH share exactly the same coalitions")


def aggregate(records,cfg,catalog):
    records=sorted(records,key=lambda r:r["scenario"]["seed"])
    if len(records)!=cfg["task_a"]["construction_scenarios"]:raise ValueError("Incomplete construction ensemble")
    if len({r["scenario"]["seed"] for r in records})!=len(records):raise ValueError("Duplicate construction seed")
    checkpoints=sorted(set(n for n in cfg["task_a"]["checkpoints"]+[len(records)] if n<=len(records)))
    result=dict(method="conditional scenario means; unseen assets explicitly fall back to CEN only when configured",checkpoints={})
    for n in checkpoints:
        sums=defaultdict(lambda:defaultdict(float));counts=defaultdict(lambda:defaultdict(int))
        for r in records[:n]:
            for strategy,scores in r["scores"].items():
                for asset,value in scores.items():sums[strategy][asset]+=value;counts[strategy][asset]+=1
        tables={s:{a:v/counts[s][a] for a,v in values.items()} for s,values in sums.items()}
        tables["CEN"]=catalog["centrality"]
        coverage={s:dict(scored_assets=len(tables[s]),missing_scores=len(catalog["assets"])-len(tables[s]),sampled_assets=len(counts[s]),total_assets=len(catalog["assets"]),unseen_assets=len(catalog["assets"])-len(counts[s]),
            minimum_observed_count=min(counts[s].values(),default=0)) for s in tables}
        result["checkpoints"][str(n)]=dict(tables=tables,counts={s:dict(v) for s,v in counts.items()},coverage=coverage)
    return result


def order(assets,strategy,tables,catalog,cfg,od_scores=None):
    power_rule=road_rule=strategy; missing=[]
    if strategy.startswith("OD_"):power_rule=strategy[3:];road_rule="OD"
    if strategy=="MIX_IJSHroad_JSHpower":power_rule="JSH";road_rule="IJSH"
    if strategy=="MIX_JSHroad_IJSHpower":power_rule="IJSH";road_rule="JSH"
    scores={}
    for a in assets:
        rule=power_rule if a.startswith("power:") else road_rule
        if rule=="OD":scores[a]=(od_scores or {}).get(a,0.0)
        elif a in tables[rule]:scores[a]=tables[rule][a]
        else:
            if cfg["task_a"]["unseen_asset_policy"]!="centrality_fallback_recorded":raise ValueError(f"Missing {rule} score for {a}")
            missing.append(a);scores[a]=catalog["centrality"][a]
    return sorted(assets,key=lambda a:(-scores[a],a)),missing


def critical_set(catalog,tag):
    powers=sorted(k for k,a in catalog["assets"].items() if a["kind"]=="power")
    critical=list(catalog["critical_substations"]); destinations=[catalog["shelter"]]
    pair_mode="all_pairs" if tag=="allpairs" else "depot_to_critical"
    if tag=="central":critical=[a.split(":",1)[1] for a in sorted(powers,key=lambda a:(-catalog["centrality"][a],a))[:2]]
    elif tag=="roadcentral":
        scores=defaultdict(float)
        for a,info in catalog["assets"].items():
            if info["kind"]=="road":
                for node in info["targets"]:scores[node]+=catalog["centrality"][a]
        destinations=sorted(scores,key=lambda n:(-scores[n],n))[:3];critical=[]
    elif tag.startswith("rand") or tag.startswith("size"):
        n=int(tag[4:]) if tag.startswith("size") else 2
        seed=20260910+n if tag.startswith("size") else {"randA":20260911,"randB":20260912,"randC":20260913}[tag]
        critical=[a.split(":",1)[1] for a in random.Random(seed).sample(powers,n)]
    elif tag=="allnodes":
        critical=[];destinations=sorted({n for a in catalog["assets"].values() if a["kind"]=="road" for n in a["targets"]})
    destinations += [catalog["assets"]["power:"+sid]["targets"][0] for sid in critical]
    destinations=sorted(set(destinations)-{catalog["depot"]})
    return dict(critical_substations=critical,destinations=destinations,pair_mode=pair_mode,tag=tag)


def od_table(catalog,cfg,fingerprint,tag="default",k=None):
    k=k or cfg["task_c"]["k"]; selection=critical_set(catalog,tag)
    key=digest(dict(selection=selection,k=k,weights=cfg["task_c"]["path_weights"]))
    path=local_path(cfg,"cache")/fingerprint/"od"/(key+".json")
    def compute():
        # Splitting each arc into its own intermediate node preserves parallel arcs
        # for Yen's loopless-path algorithm, which does not support MultiDiGraph.
        graph=nx.DiGraph(); owner={}
        for asset,info in catalog["assets"].items():
            if info["kind"]=="road":owner.update({i:asset for i in info["link_ids"]})
        for r in rows(AUSTIN/"data/processed/road/links.csv"):
            if r["centroid_connector"]=="1":continue
            u,v,i=int(r["from_node"]),int(r["to_node"]),int(r["link_id"]);cost=float(r["free_flow_time_source"])
            midpoint=("arc",i);graph.add_edge(u,midpoint,weight=cost/2);graph.add_edge(midpoint,v,weight=cost/2)
        depot=catalog["depot"];destinations=selection["destinations"]
        pairs=([(o,d) for o in [depot]+destinations for d in [depot]+destinations if o!=d]
            if selection["pair_mode"]=="all_pairs" else [(depot,d) for d in destinations])
        weights=cfg["task_c"]["path_weights"];scores=defaultdict(float);unreachable=[];path_count=0
        for origin,destination in pairs:
            try:
                paths=itertools.islice(nx.shortest_simple_paths(graph,origin,destination,weight="weight"),k)
                for rank,path_nodes in enumerate(paths):
                    weight=weights[min(rank,len(weights)-1)];path_count+=1
                    for node in path_nodes:
                        if isinstance(node,tuple) and node[1] in owner:scores[owner[node[1]]]+=weight
            except (nx.NetworkXNoPath,nx.NodeNotFound):unreachable.append([origin,destination])
        return dict(scores=dict(scores),selection=selection,k=k,od_pairs=len(pairs),paths=path_count,unreachable_pairs=unreachable,
            score_method="sum path-rank weights across arc identities, then aggregate to physical road groups")
    return cached_json(path,compute)
