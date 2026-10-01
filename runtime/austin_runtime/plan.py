"""Declarative replication matrix. Building a plan never invokes a solver."""
from __future__ import annotations
from .common import digest

STAGES=["construct","tables","a-main","a-shift","a-alpha","a-closure","a-greedy","a-cases",
        "b-main","b-sensitivity","c-main","c-k","c-decomp","c-shift","c-alternatives","analyze"]
A_RULES=["CEN","JSH","IJSH"]
C_RULES=A_RULES+["OD_CEN","OD_JSH","OD_IJSH"]
B_RULES=["baseline_reference","single_triangle","single_gini_restore","single_maximin_time_avg_cri",
         "single_p90_access_restore","weighted_gini_restore_l050","weighted_maximin_time_avg_cri_l050","guardrail_gini_restore"]
SHIFTS=["small","large","light","clustered"]


def job(stage,index,strategy,variant="main",tag="default",k=None,options=None):
    value=dict(stage=stage,index=index,strategy=strategy,variant=variant,tag=tag,k=k,options=options or {})
    value["id"]=f"{variant}_{index+1:05d}_{strategy}_{digest(value)[:12]}"
    return value


def jobs(cfg,stage,cases=None):
    a,b,c=cfg["task_a"],cfg["task_b"],cfg["task_c"]
    if stage=="construct":return [job(stage,i,"construct") for i in range(a["construction_scenarios"])]
    if stage in ("tables","analyze"):return []
    out=[]
    if stage=="a-main":
        out=[job(stage,i,s) for i in range(a["evaluation_scenarios"]) for s in A_RULES]
    elif stage=="a-shift":
        out=[job(stage,i,s,v) for v in SHIFTS for i in range(a["shift_scenarios"]) for s in A_RULES]
    elif stage=="a-alpha":
        out=[job(stage,i,"IJSH" if x==1 else f"IJSH_a{x:g}") for x in a["alpha_values"] for i in range(a["evaluation_scenarios"])]
    elif stage=="a-closure":
        out=[job(stage,i,s,options={"penalty":1000.}) for i in range(a["evaluation_scenarios"]) for s in A_RULES]
    elif stage=="a-greedy":out=[job(stage,i,s) for i in range(a["evaluation_scenarios"]) for s in A_RULES+["GREEDY"]]
    elif stage=="a-cases":
        out=[job(stage,i,s,options={"penalty":p,"case_label":label}) for label,i in (cases or [("pending-selection",0)]) for p in a["closure_penalties"] for s in A_RULES]
    elif stage=="b-main":out=[job(stage,i,s) for i in range(b["scenarios"]) for s in B_RULES]
    elif stage=="b-sensitivity":
        i=b["representative_scenario"]-1
        for s in B_RULES:
            for we,wa in b["weight_sensitivity"]:
                out.append(job(stage,i,s,options={"reevaluate":True,"electric_weight":we,"access_weight":wa}))
            for threshold in b["threshold_sensitivity"]:
                out.append(job(stage,i,s,options={"reevaluate":True,"cri_threshold":threshold,"access_threshold":threshold}))
        for metric in ("gini_restore","p90_restore","maximin_time_avg_cri_loss","p90_access_restore"):
            for lam in b["lambda_sensitivity"]:
                out.append(job(stage,i,"weighted_extra",options={"metric":metric,"lambda":lam}))
    elif stage=="c-main":out=[job(stage,i,s) for i in range(a["evaluation_scenarios"]) for s in C_RULES]
    elif stage=="c-k":
        out=[job(stage,i,s,k=k) for k in c["k_sensitivity"] for i in range(c["sensitivity_scenarios"]) for s in ("CEN","IJSH","OD_CEN","OD_IJSH")]
    elif stage=="c-decomp":
        out=[job(stage,i,s) for i in range(a["evaluation_scenarios"]) for s in A_RULES+["MIX_IJSHroad_JSHpower","MIX_JSHroad_IJSHpower"]]
    elif stage=="c-shift":
        out=[job(stage,i,s,v) for v in SHIFTS for i in range(a["shift_scenarios"]) for s in C_RULES]
    elif stage=="c-alternatives":
        for tag in c["alternative_sets"]:
            out += [job(stage,i,s,tag=tag) for i in range(a["evaluation_scenarios"]) for s in C_RULES]
            if tag in ("allpairs","central","roadcentral"):
                out += [job(stage,i,s,v,tag=tag) for v in SHIFTS for i in range(a["shift_scenarios"]) for s in C_RULES]
            if tag=="allpairs":
                out += [job(stage,i,s,tag=tag,k=k) for k in c["k_sensitivity"] for i in range(c["sensitivity_scenarios"]) for s in ("CEN","IJSH","OD_CEN","OD_IJSH")]
    else:raise ValueError(stage)
    return out
