from __future__ import annotations
import math, random
from .common import atomic_json, read_json

B_STRATEGIES=["baseline_reference","single_triangle","single_gini_restore","single_maximin_time_avg_cri",
    "single_p90_access_restore","weighted_gini_restore_l050","weighted_maximin_time_avg_cri_l050","guardrail_gini_restore"]


def base_sequence(scenario):
    sequence=sorted(scenario["damage"])
    random.Random(scenario["seed"]).shuffle(sequence)
    return sequence


def objective(metrics,spec,reference,cfg):
    kind=spec["kind"];metric=spec.get("metric","triangle_area")
    if kind=="single":return metrics[metric]
    if kind=="weighted":
        def normalized(key):
            return metrics[key]/reference[key] if reference[key]!=0 else metrics[key]
        return normalized("triangle_area")+spec["lambda"]*normalized(metric)
    if kind=="guardrail":
        limit=max(0.,reference[metric]*cfg["task_b"]["guardrail_improvement"])
        return metrics["triangle_area"]+cfg["task_b"]["guardrail_penalty"]*max(0.,metrics[metric]-limit)
    raise ValueError(kind)


def specification(strategy,options=None):
    options=options or {}
    names={"gini_restore":"gini_restore","maximin_time_avg_cri":"maximin_time_avg_cri_loss",
           "p90_access_restore":"p90_access_restore","p90_restore":"p90_restore","triangle":"triangle_area"}
    if "metric" in options:return dict(kind="weighted",metric=options["metric"],**{"lambda":options["lambda"]})
    if strategy.startswith("single_"):return dict(kind="single",metric=names[strategy[7:]])
    if strategy.startswith("weighted_"):
        name,_,lam=strategy[9:].rpartition("_l")
        return dict(kind="weighted",metric=names[name],**{"lambda":int(lam)/100})
    if strategy.startswith("guardrail_"):return dict(kind="guardrail",metric=names[strategy[10:]])
    return dict(kind="single",metric="triangle_area")


def tuples(value):
    return tuple(tuples(x) for x in value) if isinstance(value,list) else value


def anneal(engine,scenario,cfg,strategy,checkpoint,signature,options=None):
    baseline=engine.simulate(scenario,base_sequence(scenario))
    if strategy=="baseline_reference":return baseline,dict(kind="baseline",optimized=False)
    spec=specification(strategy,options);reference=baseline["metrics"]
    rng=random.Random(cfg["task_b"]["sa_seed"]);count=cfg["task_b"]["sa_iterations"]
    state=dict(signature=signature,iteration=0,current=baseline["sequence"],current_metrics=reference,
        best=baseline["sequence"],best_metrics=reference,trace=[],rng_state=rng.getstate())
    if checkpoint.exists():
        state=read_json(checkpoint)
        if state["signature"]!=signature:raise ValueError("SA checkpoint fingerprint mismatch")
        rng.setstate(tuples(state["rng_state"]))
    for i in range(state["iteration"],count):
        trial=state["current"].copy();a,b=rng.sample(range(len(trial)),2);trial[a],trial[b]=trial[b],trial[a]
        result=engine.simulate(scenario,trial);cost=objective(result["metrics"],spec,reference,cfg)
        current_cost=objective(state["current_metrics"],spec,reference,cfg)
        best_cost=objective(state["best_metrics"],spec,reference,cfg)
        temperature=cfg["task_b"]["temperature"]*cfg["task_b"]["cooling"]**i
        accept=cost<=current_cost or rng.random()<math.exp(min(0.,(current_cost-cost)/max(temperature,1e-12)))
        if accept:state.update(current=trial,current_metrics=result["metrics"])
        if cost<best_cost:state.update(best=trial,best_metrics=result["metrics"])
        state["trace"].append(dict(iteration=i+1,objective=cost,accepted=accept,best=min(cost,best_cost),temperature=temperature))
        state.update(iteration=i+1,rng_state=rng.getstate());atomic_json(checkpoint,state)
    result=engine.simulate(scenario,state["best"])
    audit=dict(specification=spec,reference=reference,iterations=count,trace=state["trace"],best_objective=objective(result["metrics"],spec,reference,cfg))
    if spec["kind"]=="guardrail":
        limit=max(0.,reference[spec["metric"]]*cfg["task_b"]["guardrail_improvement"])
        audit.update(guardrail_limit=limit,guardrail_violation=max(0.,result["metrics"][spec["metric"]]-limit))
    return result,audit
