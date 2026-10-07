from __future__ import annotations
import importlib.util, math
from collections import Counter
from pathlib import Path
from .common import AUSTIN, RUNTIME, atomic_json, cpu_budget, digest, fingerprint, local_path, read_json, sha_file
from .progress import phase


def verify_inputs(cfg,catalog):
    expected=digest(dict(traffic=cfg["traffic"],recovery=cfg["recovery"]))
    if catalog["provenance"]["prepared_config"]!=expected:
        raise ValueError("Preparation settings changed; prepare this configuration first")
    for name,expected in catalog["input_hashes"].items():
        path=local_path(cfg,"prepared")/"load_catalog.csv.gz" if name=="runtime/prepared/load_catalog.csv.gz" else AUSTIN/name
        if not path.is_file() or sha_file(path)!=expected:raise ValueError(f"Prepared input changed or missing: {path}; prepare again into a new experiment")
    for name,expected in catalog["derived_model_hashes"].items():
        path=local_path(cfg,"prepared")/name
        if not path.is_file() or sha_file(path)!=expected:raise ValueError(f"Prepared circuit changed: {path}")
    if not (RUNTIME/"build/tap-b/bin/tap").is_file():raise FileNotFoundError("Native TAP-B build is missing; run setup.sh on this host")
    return True


def doctor(cfg):
    missing=[name for name in ("numpy","scipy","networkx","opendssdirect","matplotlib") if importlib.util.find_spec(name) is None]
    catalog=local_path(cfg,"prepared")/"catalog.json"
    report=dict(dependencies_missing=missing,prepared_catalog_exists=catalog.exists(),
        native_tapb_exists=(RUNTIME/"build/tap-b/bin/tap").exists(),solvers_executed=False)
    if catalog.exists():
        try:report["input_checksums_valid"]=verify_inputs(cfg,read_json(catalog))
        except (ValueError,FileNotFoundError) as exc:report.update(input_checksums_valid=False,error=str(exc))
    return report


def validate(cfg,fp=None,budget=None):
    from .simulation import Coupled
    from .traffic import TrafficState
    import numpy as np
    catalog=read_json(local_path(cfg,"prepared")/"catalog.json");verify_inputs(cfg,catalog)
    fp=fp or fingerprint(cfg,catalog);output=local_path(cfg,"output")
    budget=budget or cpu_budget(cfg)
    report=dict(fingerprint=fp,passed=False,status="running",checks=[],
        regional_processes=budget["validation_power_workers"],
        scope="destination-host runtime acceptance: physical AC limits, coupling, TAP-B, full-damage and repair round trip")
    engine=None
    def mark(label):
        report["phase"]=label;atomic_json(output/"validation.json",report)
    try:
        mark("healthy_baseline")
        with phase(f"validation healthy baseline; AC processes={budget['validation_power_workers']}"):
            engine=Coupled(cfg,catalog,fp,output/"validation_scratch",power_workers=budget["validation_power_workers"],verbose=True)
            healthy=engine.state({})
        if not math.isclose(healthy["power_func"],1.) or not math.isclose(healthy["road_func"],1.):raise ValueError("Healthy normalization failed")
        if len(healthy["power"]["regions"])!=6 or not all(r["feasible"] and r["served_kw"]>0 for r in healthy["power"]["regions"]):
            raise ValueError("Every one of the six healthy regions must supply positive load within AC limits")
        report["checks"].append("all six original regional circuits AC-feasible with explicitly recorded curtailment")
        healthy_fraction=healthy["power"]["served_kw"]/healthy["power"]["nominal_kw"]
        if cfg["power"].get("operator","uniform_grid")=="local":
            # Limits referenced to the healthy base case: the undamaged grid must serve its full nominal load.
            if healthy_fraction<0.999:raise ValueError(f"Local operator sheds the healthy grid: served fraction {healthy_fraction:.6f}")
            report["checks"].append("healthy base case serves nominal load under base-case-referenced limits")
        counts=Counter(r["substation_id"] for r in catalog["signals"])
        sub=sorted(counts,key=lambda s:(-counts[s],s))[0];asset="power:"+sub
        mark("full_substation_fault")
        with phase(f"validation full fault: {asset}"):
            damaged=engine.state({asset:0.0})
        if any(not r["feasible"] for r in damaged["power"]["regions"]):raise ValueError("Infeasible damaged state")
        region=catalog["assets"][asset]["region"]
        # Force a new native compile/solve for the restoration round trip (bypass cache).
        mark("fresh_restoration")
        with phase(f"validation fresh restoration: {region}"):
            restored=engine.power.fresh_region(region,{})
        reference=next(r for r in healthy["power"]["regions"] if r["region"]==region)
        if not math.isclose(restored["served_kw"],reference["served_kw"],rel_tol=1e-6,abs_tol=.1):raise ValueError("Restoration depends on previous electrical state")
        if np.array_equal(healthy["traffic"].factors,damaged["traffic"].factors):
            report["checks"].append("selected full substation fault caused no signal-factor change; inspect native tie restoration/baseline shedding")
        else:report["checks"].append("AC-derived signal service changed TAP-B capacities")
        mark("partial_derating")
        with phase(f"validation partial derating: {asset}"):
            partial=engine.power.evaluate({asset:.5})
        if not all(r["feasible"] for r in partial["regions"]):raise ValueError("Partial transformer derating is infeasible")
        if cfg["power"].get("operator","uniform_grid")=="local":
            # Damage may never raise supply above the undamaged grid (the legacy uniform grid could).
            for label,served in (("full fault",damaged["power"]["served_kw"]),("partial derating",partial["served_kw"])):
                if served>healthy["power"]["served_kw"]*(1+1e-9)+1e-6:raise ValueError(f"Non-monotone power operator: {label} serves more than the healthy grid")
            report["checks"].append("full fault and partial derating serve no more than the healthy grid")
        report.update(passed=True,status="passed",phase="complete",healthy_served_kw=healthy["power"]["served_kw"],nominal_kw=healthy["power"]["nominal_kw"],
            healthy_served_fraction=healthy_fraction,power_operator=cfg["power"].get("operator","uniform_grid"),
            healthy_power_regions=healthy["power"]["regions"],healthy_traffic=healthy["traffic"].report,
            tested_substation=sub,full_fault_served_kw=damaged["power"]["served_kw"],partial_fault_served_kw=partial["served_kw"],
            damaged_traffic=damaged["traffic"].report,restoration_round_trip_served_kw=restored["served_kw"],
            modeled_signals=len(catalog["signals"]),modeled_zones=len(engine.zones),excluded_zones=engine.excluded_zones)
        report["checks"] += ["native full-fault isolation and partial thermal derating solved", "fresh physical restoration agrees with healthy state", "TAP-B row identities, relative gap and nodal flow conservation passed"]
        return report
    except BaseException as exc:
        report.update(status="failed",error=f"{type(exc).__name__}: {exc}")
        raise
    finally:
        try:atomic_json(output/"validation.json",report)
        finally:
            if engine is not None:engine.close(cancel=not report["passed"])
