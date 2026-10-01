from __future__ import annotations
import concurrent.futures as futures
import importlib.metadata, json, multiprocessing, os, resource, signal, sys, time, traceback
from collections import defaultdict
from pathlib import Path
from .common import RUNTIME, atomic_json, cpu_budget, digest, file_lock, fingerprint, local_path, read_json, sha_file
from .plan import STAGES, jobs


def result_path(output,job):return Path(output)/"results"/job["stage"]/(job["id"]+".json")


def table_required(job):return job["stage"] not in ("construct","b-main","b-sensitivity") and job["strategy"]!="GREEDY"


def signature(fp,job,output):
    table=Path(output)/"tables.json"
    if table_required(job) and not table.exists():raise FileNotFoundError("Run construct and tables stages first")
    return digest(dict(run=fp,job=job,table_sha256=sha_file(table) if table_required(job) else None))


def worker_init():
    # A process group owns its TAP-B descendants so cancellation cannot leave solvers behind.
    if hasattr(os,"setsid"):os.setsid()
    for name in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","NUMEXPR_NUM_THREADS"):
        os.environ[name]="1"


def execute(cfg,fp,job,expected_signature):
    from .simulation import Coupled
    from .scenarios import generate
    from .rankings import construct, order, od_table, critical_set
    from .experiments import anneal, base_sequence
    output=local_path(cfg,"output");catalog=read_json(local_path(cfg,"prepared")/"catalog.json")
    path=result_path(output,job)
    with file_lock(path.with_suffix(".lock")):
        if path.exists():
            saved=read_json(path)
            if saved["signature"]!=expected_signature:raise ValueError("Existing result has different inputs")
            return dict(id=job["id"],status="resumed")
        scratch=output/"scratch"/job["stage"]/job["id"];scratch.mkdir(parents=True,exist_ok=True)
        checkpoint=output/"checkpoints"/job["stage"]/(job["id"]+".json")
        previous=Path.cwd();started=time.monotonic()
        try:
            os.chdir(scratch)
            seed0=(cfg["task_a"]["construction_seed"] if job["stage"]=="construct" else
                   cfg["task_b"]["seed"] if job["stage"].startswith("b-") else
                   cfg["task_a"]["evaluation_seed"] if job["variant"]=="main" else None)
            scenario=generate(catalog,cfg,job["variant"],job["index"],seed0)
            engine=Coupled(cfg,catalog,fp,scratch)
            if job["stage"]=="construct":
                payload=construct(engine,scenario,cfg,checkpoint,expected_signature)
            elif job["stage"].startswith("b-"):
                if job["options"].get("reevaluate"):
                    refs=[r for r in results(output,"b-main",fp) if r["job"]["index"]==job["index"] and r["job"]["strategy"]==job["strategy"]]
                    if len(refs)!=1:raise ValueError("Representative Task B nominal result missing or duplicated")
                    eq={**cfg["equity"],**{k:v for k,v in job["options"].items() if k!="reevaluate"}}
                    run=engine.simulate(scenario,refs[0]["payload"]["sequence"],equity=eq)
                    audit=dict(reevaluated_fixed_sequence=True,equity=eq)
                else:run,audit=anneal(engine,scenario,cfg,job["strategy"],checkpoint,expected_signature,job["options"])
                payload={**run,"optimization":audit,"fallback_assets":[]}
            else:
                table_data=read_json(output/"tables.json") if table_required(job) else None
                tables=table_data["checkpoints"][str(max(map(int,table_data["checkpoints"])))]["tables"] if table_data else {}
                od=None;critical=catalog["critical_substations"]
                if job["stage"].startswith("c-"):
                    od=od_table(catalog,cfg,fp,job["tag"],job["k"])
                    critical=od["selection"]["critical_substations"] or catalog["critical_substations"]
                if job["strategy"]=="GREEDY":sequence=sorted(scenario["damage"]);missing=[]
                else:sequence,missing=order(scenario["damage"],job["strategy"],tables,catalog,cfg,od["scores"] if od else None)
                run=engine.simulate(scenario,sequence,greedy=job["strategy"]=="GREEDY",penalty=job["options"].get("penalty"),critical=critical)
                payload={**run,"fallback_assets":missing,"od_metadata":{k:v for k,v in od.items() if k!="scores"} if od else None,
                    "od_damaged_roads":sum(a.startswith("road:") for a in scenario["damage"]) if od else None,
                    "od_damaged_roads_nonzero":sum(a.startswith("road:") and od["scores"].get(a,0)>0 for a in scenario["damage"]) if od else None}
            record=dict(signature=expected_signature,fingerprint=fp,job=job,payload=payload,
                resources=dict(wall_seconds=time.monotonic()-started,worker_lifetime_peak_rss_mib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,
                    native_children_peak_rss_mib=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss/1024))
            atomic_json(path,record)
            error=path.with_suffix(".error.json")
            if error.exists():error.unlink()
            return dict(id=job["id"],status="complete")
        except Exception as exc:
            atomic_json(path.with_suffix(".error.json"),dict(signature=expected_signature,fingerprint=fp,job=job,
                exception=type(exc).__name__,message=str(exc),traceback=traceback.format_exc()))
            return dict(id=job["id"],status="failed",error=str(exc))
        finally:os.chdir(previous)


def results(output,stage,fp):
    for path in sorted((Path(output)/"results"/stage).glob("*.json")):
        if path.name.endswith(".error.json"):continue
        record=read_json(path)
        if record["fingerprint"]!=fp:raise ValueError("Mixed result fingerprints")
        yield record


def select_cases(output,cfg,fp):
    by_case=defaultdict(dict)
    for r in results(output,"a-main",fp):by_case[r["job"]["index"]][r["job"]["strategy"]]=r["payload"]["metrics"]["triangle_area"]
    if len(by_case)!=cfg["task_a"]["evaluation_scenarios"] or any(set(r)!={"CEN","JSH","IJSH"} for r in by_case.values()):
        raise ValueError("Representative cases require complete Task A main results")
    pairs=sorted(((v["JSH"]-v["IJSH"])/max(abs(v["JSH"]),1e-12),i) for i,v in by_case.items())
    cases=[("typical",pairs[len(pairs)//2][1]),("largest_gain",pairs[-1][1]),("largest_loss",pairs[0][1])]
    cases=cases[:cfg["task_a"]["representative_cases"]]
    atomic_json(Path(output)/"case_selection.json",dict(rule="median, maximum, minimum paired JSH-to-IJSH relative improvement; selected before closure reruns",cases=cases))
    return cases


def run_stage(cfg,stage,fp,budget):
    from .rankings import aggregate
    output=local_path(cfg,"output");catalog=read_json(local_path(cfg,"prepared")/"catalog.json")
    if stage=="tables":
        records=[r["payload"] for r in results(output,"construct",fp)]
        data=aggregate(records,cfg,catalog);data["fingerprint"]=fp
        atomic_json(output/"tables.json",data);return
    if stage=="analyze":
        from .analysis import analyze
        analyze(cfg,fp);return
    cases=select_cases(output,cfg,fp) if stage=="a-cases" else None
    todo=jobs(cfg,stage,cases)
    atomic_json(output/"plans"/(stage+".json"),dict(fingerprint=fp,jobs=todo))
    pending=[]
    for item in todo:
        sig=signature(fp,item,output);path=result_path(output,item)
        if path.exists():
            if read_json(path)["signature"]!=sig:raise ValueError(f"Stale result: {path}; choose a new output directory")
            stale_error=path.with_suffix(".error.json")
            if stale_error.exists():stale_error.unlink()
        else:pending.append((item,sig))
    print(f"{stage}: {len(pending)} pending / {len(todo)} total; {budget['workers']} workers",flush=True)
    if not pending:
        failures=output/"plans"/(stage+".failures.json")
        if failures.exists():failures.unlink()
        return
    errors=[];completed=0
    wave_size=budget["workers"]*cfg["runtime"]["jobs_per_worker_batch"]
    # Recreate complete pools at bounded batch boundaries. Avoid the known
    # max_tasks_per_child deadlock in older supported CPython versions.
    for offset in range(0,len(pending),wave_size):
        pool=futures.ProcessPoolExecutor(max_workers=budget["workers"],
            mp_context=multiprocessing.get_context("spawn"),initializer=worker_init)
        try:
            iterator=iter(pending[offset:offset+wave_size]);running={}
            def submit_one():
                try:item,sig=next(iterator)
                except StopIteration:return False
                running[pool.submit(execute,cfg,fp,item,sig)]=item;return True
            for _ in range(budget["workers"]):submit_one()
            while running:
                done,_=futures.wait(running,return_when=futures.FIRST_COMPLETED)
                for future in done:
                    item=running.pop(future)
                    try:status=future.result()
                    except Exception as exc:
                        status=dict(id=item["id"],status="failed",error=str(exc))
                        atomic_json(result_path(output,item).with_suffix(".error.json"),
                            dict(fingerprint=fp,job=item,message=str(exc),exception=type(exc).__name__,traceback=traceback.format_exc()))
                    completed+=1;print(f"{stage} {completed}/{len(pending)}: {status}",flush=True)
                    if status["status"]=="failed":errors.append(status)
                    submit_one()
            pool.shutdown(wait=True)
        except BaseException:
            # Linux workers own their sessions, including native TAP-B children.
            for process in list((getattr(pool,"_processes",None) or {}).values()):
                try:os.killpg(process.pid,signal.SIGTERM)
                except (ProcessLookupError,PermissionError):pass
            pool.shutdown(wait=False,cancel_futures=True)
            if errors:atomic_json(output/"plans"/(stage+".failures.json"),errors)
            raise
    if errors:
        atomic_json(output/"plans"/(stage+".failures.json"),errors)
        raise RuntimeError(f"{stage}: {len(errors)} jobs failed. Successful jobs/checkpoints are retained; rerun after resolving errors.")
    failures=output/"plans"/(stage+".failures.json")
    if failures.exists():failures.unlink()


def run(cfg,stage="all",requested=None):
    from .validation import validate, verify_inputs
    catalog_path=local_path(cfg,"prepared")/"catalog.json"
    if not catalog_path.exists():raise FileNotFoundError("Run the explicit prepare step first")
    catalog=read_json(catalog_path);verify_inputs(cfg,catalog)
    fp=fingerprint(cfg,catalog);output=local_path(cfg,"output")
    output.mkdir(parents=True,exist_ok=True);budget=cpu_budget(cfg,requested) if stage!="analyze" else {}
    with file_lock(output/"run.lock",blocking=False):
        manifest=output/"run_manifest.json"
        if manifest.exists() and read_json(manifest)["fingerprint"]!=fp:
            raise ValueError("Configuration, code or data changed; use a new output directory. Existing results were not overwritten.")
        versions={k:importlib.metadata.version(k) for k in ("numpy","scipy","networkx","OpenDSSDirect.py")}
        if stage!="analyze" or not manifest.exists():
            atomic_json(manifest,dict(fingerprint=fp,config=cfg,budget=budget,python=sys.version,packages=versions,
                source="Austin independent runtime; legacy results never imported",validation_required=True))
        if stage!="analyze":
            gate=output/"validation.json"
            if not gate.exists() or read_json(gate).get("fingerprint")!=fp or not read_json(gate).get("passed"):
                validate(cfg,fp)
        for selected in (STAGES if stage=="all" else [stage]):run_stage(cfg,selected,fp,budget)
    return dict(output=str(output),fingerprint=fp,stage=stage,status="completed")
