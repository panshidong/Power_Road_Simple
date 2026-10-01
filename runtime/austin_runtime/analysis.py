from __future__ import annotations
import math
from collections import defaultdict
from pathlib import Path
import numpy as np
from scipy.stats import binomtest, kendalltau, spearmanr
from .common import atomic_json, digest, local_path, read_json, write_rows


def finite(value):return float(value) if math.isfinite(float(value)) else None


def comparison(reference,candidate,cfg,seed,higher_is_better=False):
    a=np.asarray(reference,dtype=float);b=np.asarray(candidate,dtype=float);delta=(b-a) if higher_is_better else (a-b)
    tol=1e-8;wins=int(sum(delta>tol));losses=int(sum(delta<-tol));ties=len(delta)-wins-losses
    rng=np.random.default_rng(seed);draws=cfg["analysis"]["bootstrap_resamples"];trim=cfg["analysis"]["trim_each_tail"]
    boot=[];trimmed=[]
    for start in range(0,draws,256):
        sample=delta[rng.integers(0,len(delta),size=(min(256,draws-start),len(delta)))]
        boot.extend(np.mean(sample,axis=1).tolist());sample.sort(axis=1);cut=int(sample.shape[1]*trim)
        trimmed.extend(np.mean(sample[:,cut:-cut] if cut else sample,axis=1).tolist())
    sorted_delta=np.sort(delta);cut=int(len(delta)*trim)
    return dict(n=len(delta),higher_is_better=higher_is_better,positive_difference="candidate improvement",reference_mean=float(a.mean()),candidate_mean=float(b.mean()),
        mean_reduction=float(delta.mean()),median_reduction=float(np.median(delta)),
        percent_reduction_of_reference_mean=float(100*delta.mean()/a.mean()) if a.mean()!=0 else None,
        mean_reduction_ci95=np.quantile(boot,[.025,.975]).tolist(),
        trimmed_mean_reduction=float(np.mean(sorted_delta[cut:-cut] if cut else sorted_delta)),
        trimmed_mean_ci95=np.quantile(trimmed,[.025,.975]).tolist(),wins=wins,losses=losses,ties=ties,
        exact_sign_p=float(binomtest(wins,wins+losses,.5).pvalue) if wins+losses else 1.0)


def analyze(cfg,fp):
    output=local_path(cfg,"output");dest=output/"analysis";dest.mkdir(parents=True,exist_ok=True)
    flat=[];groups=defaultdict(lambda:defaultdict(dict)); failures=[];expected=0;missing=[]
    for plan in sorted((output/"plans").glob("*.json")):
        if plan.name.endswith(".failures.json"):continue
        manifest=read_json(plan)
        if manifest.get("fingerprint")!=fp:raise ValueError("Mixed plan fingerprints in analysis")
        for job in manifest.get("jobs",[]):
            if job["stage"]=="construct":continue
            expected+=1
            if not (output/"results"/job["stage"]/(job["id"]+".json")).exists():missing.append(job)
    for path in sorted((output/"results").glob("*/*.json")):
        record=read_json(path)
        if record.get("fingerprint")!=fp:raise ValueError("Mixed result fingerprints in analysis")
        if path.name.endswith(".error.json"):failures.append(dict(path=str(path),message=record["message"]));continue
        job=record["job"];run=record["payload"]
        if job["stage"]=="construct":continue
        group=(job["stage"],job["variant"],job["tag"],job["k"],digest(job["options"])[:12])
        # Keep aggregate inputs only; full event arrays stay on disk even for large ensembles.
        groups[group][job["strategy"]][run["scenario_seed"]]={k:run[k] for k in ("metrics","complete","fallback_assets","sequence")}

        flat.append(dict(stage=job["stage"],variant=job["variant"],tag=job["tag"],k=job["k"],options_digest=group[-1],
            strategy=job["strategy"],scenario_id=run["scenario_id"],scenario_seed=run["scenario_seed"],complete=run["complete"],
            fallback_assets=len(run.get("fallback_assets",[])),healthy_served_fraction=run["healthy_served_fraction"],
            od_damaged_roads=run.get("od_damaged_roads"),od_damaged_roads_nonzero=run.get("od_damaged_roads_nonzero"),
            case_label=job["options"].get("case_label",""),closure_penalty=job["options"].get("penalty",cfg["traffic"]["closure_penalty"]),
            **record.get("resources",{}),**run["metrics"]))
    write_rows(dest/"scenario_rows.csv",flat)
    paired=[];summaries=[];agreement=[]
    for group,strategies in groups.items():
        for strategy,records in strategies.items():
            for metric in ("triangle_area","weighted_triangle_area","gini_restore","min_time_avg_cri","p90_access_restore"):
                values=[r["metrics"][metric] for r in records.values()]
                summaries.append(dict(stage=group[0],variant=group[1],tag=group[2],k=group[3],options_digest=group[4],strategy=strategy,
                    metric=metric,n=len(values),mean=float(np.mean(values)),median=float(np.median(values)),p95=float(np.quantile(values,.95)),
                    incomplete_recoveries=sum(not r["complete"] for r in records.values()),
                    runs_using_unseen_score_fallback=sum(bool(r.get("fallback_assets")) for r in records.values())))
        pairs=[("CEN","JSH"),("JSH","IJSH"),("CEN","IJSH"),("JSH","MIX_IJSHroad_JSHpower"),("JSH","MIX_JSHroad_IJSHpower")]
        pairs += [(s,"OD_"+s) for s in ("CEN","JSH","IJSH")]
        pairs += [(s,"OD_CEN") for s in ("JSH","IJSH")]+[(s,"GREEDY") for s in ("CEN","JSH","IJSH")]
        pairs += [("baseline_reference",s) for s in strategies if s!="baseline_reference"]
        pairs += [("IJSH",s) for s in strategies if s.startswith("IJSH_a")]
        for ref,candidate in pairs:
            if ref not in strategies or candidate not in strategies:continue
            matched=sorted(set(strategies[ref])&set(strategies[candidate]))
            if not matched:continue
            for kind in ("power:","road:"):
                same=sum([a for a in strategies[ref][seed]["sequence"] if a.startswith(kind)]==[a for a in strategies[candidate][seed]["sequence"] if a.startswith(kind)] for seed in matched)
                agreement.append(dict(stage=group[0],variant=group[1],tag=group[2],k=group[3],options_digest=group[4],reference=ref,candidate=candidate,asset_type=kind,n=len(matched),identical_sequences=same,identical_fraction=same/len(matched)))
                if candidate=="OD_"+ref and kind=="power:" and same!=len(matched):
                    raise ValueError("OD paired comparison changed the power ranking")
            for metric in ("triangle_area","weighted_triangle_area","gini_restore","p90_access_restore","min_time_avg_cri"):
                a=[strategies[ref][seed]["metrics"][metric] for seed in matched]
                b=[strategies[candidate][seed]["metrics"][metric] for seed in matched]
                seed=int(digest([cfg["analysis"]["seed"],group,ref,candidate,metric])[:16],16)
                stats=comparison(a,b,cfg,seed,higher_is_better=metric=="min_time_avg_cri")
                paired.append(dict(stage=group[0],variant=group[1],tag=group[2],k=group[3],options_digest=group[4],reference=ref,candidate=candidate,metric=metric,
                    reference_unpaired=len(strategies[ref])-len(matched),candidate_unpaired=len(strategies[candidate])-len(matched),**stats))
    write_rows(dest/"aggregate.csv",summaries);atomic_json(dest/"paired_statistics.json",paired)
    write_rows(dest/"sequence_agreement.csv",agreement)
    od_scores=[];od_coverage=[]
    for path in sorted((local_path(cfg,"cache")/fp/"od").glob("*.json")):
        table=read_json(path);selection=table["selection"]
        od_coverage.append(dict(tag=selection["tag"],k=table["k"],od_pairs=table["od_pairs"],paths=table["paths"],unreachable_pairs=len(table["unreachable_pairs"]),positive_road_groups=sum(v>0 for v in table["scores"].values())))
        od_scores.extend(dict(tag=selection["tag"],k=table["k"],asset=asset,score=score) for asset,score in table["scores"].items())
    write_rows(dest/"od_coverage.csv",od_coverage);write_rows(dest/"od_asset_scores.csv",od_scores)
    stability=[];table_file=output/"tables.json"
    if table_file.exists():
        tables=read_json(table_file);last=tables["checkpoints"][str(max(map(int,tables["checkpoints"])))]
        for n,checkpoint in tables["checkpoints"].items():
            for strategy,values in checkpoint["tables"].items():
                for kind in ("power:","road:"):
                    shared=sorted(a for a in values if a.startswith(kind) and a in last["tables"][strategy])
                    x=[values[a] for a in shared];y=[last["tables"][strategy][a] for a in shared]
                    rho=tau=None
                    if len(shared)>1 and len(set(x))>1 and len(set(y))>1:
                        rho=finite(spearmanr(x,y).statistic);tau=finite(kendalltau(x,y).statistic)
                    stability.append(dict(checkpoint=int(n),strategy=strategy,asset_type=kind,shared_assets=len(shared),spearman=rho,kendall=tau))
        atomic_json(dest/"table_coverage.json",{n:c["coverage"] for n,c in tables["checkpoints"].items()})
        write_rows(dest/"table_stability.csv",stability)
    # Scenario-level Pareto labels keep incomparable scenarios/options separate.
    pareto=[]
    for group,strategies in groups.items():
        if group[0]!="b-main":continue
        for seed in sorted({seed for rs in strategies.values() for seed in rs}):
            points={s:rs[seed]["metrics"] for s,rs in strategies.items() if seed in rs}
            for s,m in points.items():
                dominated=any(q["triangle_area"]<=m["triangle_area"] and q["gini_restore"]<=m["gini_restore"] and
                    (q["triangle_area"]<m["triangle_area"] or q["gini_restore"]<m["gini_restore"]) for t,q in points.items() if t!=s)
                pareto.append(dict(scenario_seed=seed,strategy=s,triangle_area=m["triangle_area"],gini_restore=m["gini_restore"],pareto=not dominated))
    write_rows(dest/"task_b_pareto.csv",pareto)
    status=dict(fingerprint=fp,planned_evaluation_jobs=expected,completed_rows=len(flat),missing_jobs=missing,failures=failures,
        complete=bool(expected) and not missing and not failures,completeness_scope="stages with a saved plan; see task matrix for entire study",uncensored_rows=sum(r["complete"] for r in flat),
        interpretation="paired Monte Carlo estimates conditional on recorded Austin assumptions; no comparison to legacy article numbers")
    atomic_json(dest/"status.json",status)
    plot(dest,flat,paired,stability)
    from .figures import supplemental_figures
    supplemental_figures(output,dest,flat,stability,fp)
    summary=["# Austin runtime results", "", f"Completed {len(flat)} / {expected} planned evaluation jobs; complete={status['complete']}.",
        f"Incomplete recovery trajectories: {sum(not r['complete'] for r in flat)}. Solver failures are listed separately; never replaced with zero.",
        "", "See scenario_rows.csv, paired_statistics.json, table_coverage.json, task_b_pareto.csv and figures/.",
        "Statistics include horizon-censored recovery outcomes; inspect censored zone counts and raw event logs."]
    (dest/"SUMMARY.md").write_text("\n".join(summary)+"\n")
    return status


def plot(dest,records,paired,stability):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    folder=dest/"figures";folder.mkdir(exist_ok=True)
    for stage in ("a-main","b-main","c-main","a-closure","a-greedy"):
        subset=[r for r in records if r["stage"]==stage]
        names=sorted({r["strategy"] for r in subset})
        if not names:continue
        fig,ax=plt.subplots(figsize=(max(7,len(names)*1.25),5))
        ax.boxplot([[r["triangle_area"] for r in subset if r["strategy"]==name] for name in names],tick_labels=names,showfliers=False)
        ax.set_ylabel("System resilience-loss area (assumed minutes)");ax.tick_params(axis="x",rotation=25)
        fig.tight_layout();fig.savefig(folder/(stage+".png"),dpi=180);fig.savefig(folder/(stage+".svg"));plt.close(fig)
    subset=[r for r in records if r["stage"]=="b-main"]
    if subset:
        fig,ax=plt.subplots(figsize=(8,5))
        for s in sorted({r["strategy"] for r in subset}):
            group=[r for r in subset if r["strategy"]==s]
            ax.scatter(np.mean([r["triangle_area"] for r in group]),np.mean([r["gini_restore"] for r in group]),label=s)
        ax.set_xlabel("Mean resilience-loss area");ax.set_ylabel("Mean restoration-time Gini")
        ax.legend(fontsize=7);fig.tight_layout();fig.savefig(folder/"task_b_tradeoff.png",dpi=180);fig.savefig(folder/"task_b_tradeoff.svg");plt.close(fig)
