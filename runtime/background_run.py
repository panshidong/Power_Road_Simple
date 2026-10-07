#!/usr/bin/env python3
"""Detached Austin research supervisor with frozen inputs and explicit stop/resume."""
from __future__ import annotations

import argparse
import ctypes
import fcntl
import json
import math
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import time
import traceback
from datetime import datetime
from zoneinfo import ZoneInfo


RUNTIME=Path(__file__).resolve().parent
PROJECT=RUNTIME.parent
TERMINAL={"completed","failed","stopped","disk_limit","memory_limit","stalled","time_limit"}
CORE_STAGES=["construct","tables","a-main","b-main","c-main","analyze"]


def now():return datetime.now(ZoneInfo("America/Chicago")).isoformat(timespec="seconds")


def read(path):return json.loads(Path(path).read_text())


def write(path,value):
    from austin_runtime.common import atomic_json
    atomic_json(path,value)


def lock(path):
    handle=Path(path).open("a+")
    try:fcntl.flock(handle,fcntl.LOCK_EX|fcntl.LOCK_NB)
    except BlockingIOError:
        handle.close()
        raise RuntimeError(f"Another launcher/supervisor owns {path}") from None
    return handle


def proc_info(pid):
    try:
        fields=Path(f"/proc/{pid}/stat").read_text().rsplit(")",1)[1].split()
        return int(fields[1]),fields[19]
    except (FileNotFoundError,ProcessLookupError,PermissionError):return None


def descendants(root):
    processes={int(p.name):proc_info(p.name) for p in Path("/proc").iterdir() if p.name.isdigit()}
    owned={root};found={}
    while True:
        added={pid:info[1] for pid,info in processes.items() if info and info[0] in owned and pid not in owned}
        if not added:return found
        found.update(added);owned.update(added)


def signal_owned(processes,sig):
    for pid,birth in processes.items():
        current=proc_info(pid)
        if current and current[1]==birth:
            try:os.kill(pid,sig)
            except ProcessLookupError:pass


def subreaper():
    # Linux PR_SET_CHILD_SUBREAPER. Adopt orphaned native grandchildren so a
    # crashed coordinator cannot leave its separate-session workers running.
    libc=ctypes.CDLL(None,use_errno=True)
    if libc.prctl(36,1,0,0,0):raise OSError(ctypes.get_errno(),"Cannot enable child subreaper")


def stop_children(child):
    if child is not None and child.poll() is None:
        child.send_signal(signal.SIGINT)
        try:child.wait(timeout=20)
        except subprocess.TimeoutExpired:pass
    remaining=descendants(os.getpid())
    signal_owned(remaining,signal.SIGTERM)
    if child is not None:
        try:child.wait(timeout=5)
        except subprocess.TimeoutExpired:pass
    deadline=time.monotonic()+5
    while time.monotonic()<deadline:
        remaining=descendants(os.getpid())
        if not remaining:break
        for pid in remaining:
            try:os.waitpid(pid,os.WNOHANG)
            except ChildProcessError:pass
        time.sleep(.1)
    signal_owned(descendants(os.getpid()),signal.SIGKILL)
    if child is not None:
        try:child.wait(timeout=5)
        except subprocess.TimeoutExpired:pass
    while True:
        try:
            pid,_=os.waitpid(-1,os.WNOHANG)
            if not pid:break
        except ChildProcessError:break


def progress(output,cache):
    values={"results":0,"errors":0,"cache_states":0,"construction_permutations":0,"sa_iterations":0}
    for p in (output/"results").glob("*/*.json"):
        values["errors" if p.name.endswith(".error.json") else "results"]+=1
    for p in (output/"checkpoints").glob("*/*.json"):
        data=read(p)
        values["construction_permutations"]+=data.get("permutations",0)
        values["sa_iterations"]+=data.get("iteration",0)
    for folder in [cache/"traffic",cache/"od",*(cache/"power").glob("*")]:
        if folder.is_dir():
            values["cache_states"]+=sum(p.suffix in (".json",".npz") for p in folder.iterdir())
    gate=output/"validation.json"
    if gate.exists():
        validation=read(gate)
        values["validation_phase"]=validation.get("phase")
        values["validation_passed"]=validation.get("passed",False)
    return values


def available_memory():
    for line in Path("/proc/meminfo").read_text().splitlines():
        if line.startswith("MemAvailable:"):return int(line.split()[1])*1024
    raise RuntimeError("MemAvailable unavailable")


def stage_entry(output,stage):
    from austin_runtime.common import fingerprint,read_json
    from austin_runtime.workflow import run
    spec=read(output/"control/launch.json")
    cfg=spec["config"]
    actual=fingerprint(cfg,read_json(Path(cfg["runtime"]["prepared"])/"catalog.json"))
    if actual!=spec["fingerprint"]:raise RuntimeError("Frozen source/dependencies/config fingerprint changed; refusing mixed results")
    run(cfg,stage,requested=spec["workers"])


def check_snapshot(output):
    from austin_runtime.common import fingerprint,read_json
    from austin_runtime.validation import verify_inputs
    spec=read(output/"control/launch.json");cfg=spec["config"]
    catalog=read_json(Path(cfg["runtime"]["prepared"])/"catalog.json")
    verify_inputs(cfg,catalog)
    fp=fingerprint(cfg,catalog)
    if fp!=spec["fingerprint"]:raise RuntimeError("Snapshot numerical fingerprint differs from launch plan")
    checked=dict(passed=True,fingerprint=fp,checked_at=now(),runtime=str(RUNTIME),solvers_started=0)
    write(output/"control/snapshot_check.json",checked)
    print(json.dumps(checked,indent=2),flush=True)


def prepare(args):
    from austin_runtime.common import config,cpu_budget,fingerprint,read_json,sha_file
    from austin_runtime.native_env import allocator_path
    from austin_runtime.plan import jobs,STAGES
    from austin_runtime.validation import verify_inputs
    cfg=config(args.config)
    allocator=allocator_path(cfg)
    stages=STAGES if args.profile=="full" else CORE_STAGES
    prepared=(RUNTIME/cfg["runtime"]["prepared"]).resolve()
    catalog=read_json(prepared/"catalog.json")
    verify_inputs(cfg,catalog)
    budget=cpu_budget(cfg,args.workers)
    if shutil.disk_usage(PROJECT).free<args.min_free_gib*2**30+2*2**30:
        raise RuntimeError("Insufficient free disk for snapshot plus requested reserve")
    output=args.output.resolve();output.mkdir(parents=True,exist_ok=False)
    control=output/"control";control.mkdir()
    snapshot=output/"snapshot";snap_runtime=snapshot/"runtime";snap_runtime.mkdir(parents=True)
    # Copy, never hard-link mutable source/input files. Dependencies stay in the
    # installed venv and are checked through the numerical fingerprint per stage.
    shutil.copytree(PROJECT/"data",snapshot/"data")
    shutil.copytree(prepared,snap_runtime/"prepared")
    shutil.copytree(RUNTIME/"austin_runtime",snap_runtime/"austin_runtime",ignore=shutil.ignore_patterns("__pycache__","*.pyc"))
    shutil.copytree(RUNTIME/"configs",snap_runtime/"configs")
    shutil.copy2(__file__,snap_runtime/"background_run.py")
    binary=snap_runtime/"build/tap-b/bin/tap";binary.parent.mkdir(parents=True)
    shutil.copy2(RUNTIME/"build/tap-b/bin/tap",binary)
    if allocator:
        target=snap_runtime/"build/native"/allocator.name;target.parent.mkdir(parents=True)
        shutil.copy2(allocator,target)
        license_file=allocator.parent/"jemalloc-COPYRIGHT"
        if license_file.exists():shutil.copy2(license_file,target.parent/license_file.name)
        cfg["power"]["native_allocator_path"]=str(target)
    cfg["runtime"].update(output=str(output),prepared=str(snap_runtime/"prepared"),cache=str((RUNTIME/cfg["runtime"]["cache"]).resolve()))
    fp=fingerprint(cfg,catalog)
    spec=dict(created_at=now(),source_project=str(PROJECT),interpreter=sys.executable,config=cfg,
              fingerprint=fp,workers=budget["workers"],budget=budget,stages=stages,
              profile=args.profile,omitted_stages=[s for s in STAGES if s not in stages],
              max_run_hours=args.max_run_hours,
              stage_job_counts={s:len(jobs(cfg,s,[(str(i),i) for i in range(cfg["task_a"]["representative_cases"])])) for s in stages},
              min_free_gib=args.min_free_gib,min_available_gib=4,stall_hours=args.stall_hours,poll_seconds=60,
              automatic_retries=0,scientific_status="Current disaster ensembles are uncalibrated; "+
                  (f"{cfg['recovery']['horizon_minutes']:g}-minute observation horizon retained." if cfg['recovery'].get('stop_at_horizon',True)
                   else "no observation horizon: run until all repairs complete; unreachable pending work is a reported error."),
              supervisor_sha256=sha_file(snap_runtime/"background_run.py"))
    spec["native_state_recovery"]=dict(attempts=cfg["power"].get("native_recovery_attempts",0),
        independent_confirmation_required=bool(cfg["power"].get("native_recovery_attempts",0)),
        scope="Only native signal crashes inside an isolated AC state; no experiment/stage restarts")
    if args.profile=="core":
        spec["scientific_status"] += " Core pilot only: reduced Monte Carlo/optimization budgets, no sensitivity/shift/greedy claims; inspect score coverage and fallback frequency before interpreting strategy differences."
    write(control/"launch.json",spec)
    write(control/"status.json",dict(status="prepared",prepared_at=now(),completed_stages=[],fingerprint=fp))
    (output/"RUN_SCOPE.txt").write_text(
        f"Profile: {args.profile}\nSelected stages: {', '.join(stages)}\n"
        f"Omitted stages: {', '.join(spec['omitted_stages'])}\n"
        f"Active runtime cap across resumes: {args.max_run_hours or 'unlimited'} hours\n"
        f"{spec['scientific_status']}\n",encoding="utf-8")
    write(RUNTIME/"output/latest_background.json",dict(output=str(output)))
    print(json.dumps(dict(output=str(output),fingerprint=fp,workers=budget["workers"],jobs=sum(spec["stage_job_counts"].values()),status="prepared"),indent=2),flush=True)
    return output


def launch(output):
    output=output.resolve();control=output/"control"
    with lock(control/"launch.lock"):
        # Probe supervisor lock; process identity/heartbeat is not used as a lock.
        with lock(control/"supervisor.lock"):pass
        previous=read(control/"status.json")
        if previous.get("status")=="completed":raise RuntimeError("Run is already complete")
        request=control/"stop.request"
        if request.exists():
            request.rename(control/f"stop-request-{time.time_ns()}.json")
        spec=read(control/"launch.json")
        frozen=output/"snapshot/runtime/background_run.py"
        with (control/"supervisor.log").open("a") as log:
            child=subprocess.Popen([spec["interpreter"],"-u",str(frozen),"_supervise",str(output)],
                                   stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,
                                   start_new_session=True,close_fds=True,cwd=frozen.parent)
        deadline=time.monotonic()+10
        while time.monotonic()<deadline:
            status=read(control/"status.json")
            if status.get("supervisor_pid")==child.pid:
                print(json.dumps(status,ensure_ascii=False,indent=2),flush=True)
                if status["status"] in TERMINAL and status["status"]!="completed":raise RuntimeError("Supervisor stopped at launch; see status/log")
                return
            if child.poll() is not None:raise RuntimeError(f"Supervisor exited {child.returncode}; see {control/'supervisor.log'}")
            time.sleep(.1)
        raise RuntimeError("No supervisor startup acknowledgement; inspect status before retrying")


def supervise(output):
    control=output/"control";spec=read(control/"launch.json")
    with lock(control/"supervisor.lock"):
        subreaper()
        state=read(control/"status.json")
        previous_active_seconds=state.get("active_seconds",0.)
        resumed_monotonic=time.monotonic()
        limit=spec.get("max_run_hours",0)*3600
        def time_exhausted():
            state["active_seconds"]=previous_active_seconds+time.monotonic()-resumed_monotonic
            return bool(limit and state["active_seconds"]>=limit)
        state.update(status="running",supervisor_pid=os.getpid(),started_or_resumed_at=now(),heartbeat=now())
        for key in ("error","reason","finished_or_stopped_at","last_exit_code","traceback"):
            state.pop(key,None)
        write(control/"status.json",state)
        child=None
        requested=[]
        def stop_signal(signum,frame):requested.append(signal.Signals(signum).name)
        signal.signal(signal.SIGTERM,stop_signal);signal.signal(signal.SIGINT,stop_signal)
        cache=Path(spec["config"]["runtime"]["cache"])/spec["fingerprint"]
        try:
            for stage in spec["stages"]:
                if stage in state["completed_stages"]:continue
                if time_exhausted():
                    state.update(status="time_limit",reason="Cumulative active runtime budget exhausted; saved work retained");break
                if requested or (control/"stop.request").exists():
                    state.update(status="stopped",reason="Stop requested before next stage");break
                free=shutil.disk_usage(output).free
                if free<spec["min_free_gib"]*2**30:
                    state.update(status="disk_limit",reason="Free disk below reserve");break
                log=control/f"{stage}-{time.time_ns()}.log"
                state.update(stage=stage,stage_log=str(log),stage_started_at=now(),heartbeat=now())
                write(control/"status.json",state)
                print(f"{now()} START {stage}: {log}",flush=True)
                command=[spec["interpreter"],"-u",str(Path(__file__).resolve()),"_stage",str(output),"--stage",stage]
                last_progress=time.monotonic();last=None;low_memory=0
                with log.open("w") as handle:
                    child=subprocess.Popen(command,stdin=subprocess.DEVNULL,stdout=handle,stderr=subprocess.STDOUT,
                                           cwd=RUNTIME,start_new_session=True,close_fds=True)
                    while True:
                        observed=progress(output,cache)
                        if observed!=last:last=observed;last_progress=time.monotonic()
                        free=shutil.disk_usage(output).free;ram=available_memory()
                        low_memory=low_memory+1 if ram<spec["min_available_gib"]*2**30 else 0
                        state.update(heartbeat=now(),stage_pid=child.pid,progress=observed,
                                     free_disk_gib=free/2**30,available_memory_gib=ram/2**30,
                                     seconds_without_saved_progress=time.monotonic()-last_progress)
                        code=child.poll()
                        if code is not None:
                            state["last_exit_code"]=code
                            if code:
                                state.update(status="failed",reason=f"Stage {stage} exited {code}; no automatic retry")
                            else:
                                state["completed_stages"].append(stage)
                                print(f"{now()} DONE {stage}",flush=True)
                            stop_children(child);child=None
                            write(control/"status.json",state);break
                        if requested or (control/"stop.request").exists():
                            state.update(status="stopped",reason="Stop requested")
                        elif free<spec["min_free_gib"]*2**30:
                            state.update(status="disk_limit",reason="Free disk below reserve")
                        elif low_memory>=3:
                            state.update(status="memory_limit",reason="Available RAM below 4 GiB for three consecutive polls")
                        elif time_exhausted():
                            state.update(status="time_limit",reason="Cumulative active runtime budget exhausted; saved work retained")
                        elif time.monotonic()-last_progress>spec["stall_hours"]*3600:
                            state.update(status="stalled",reason="No saved result/checkpoint/cache progress within stall limit")
                        write(control/"status.json",state)
                        if state["status"]!="running":
                            stop_children(child);child=None;break
                        # Signals and stop requests remain responsive between scans.
                        deadline=time.monotonic()+spec["poll_seconds"]
                        while time.monotonic()<deadline and child.poll() is None and not requested and not (control/"stop.request").exists() and not time_exhausted():
                            time.sleep(1)
                if state["status"]!="running":break
            else:state.update(status="completed",stage=None)
        except BaseException as exc:
            state.update(status="failed",error=f"{type(exc).__name__}: {exc}",traceback=traceback.format_exc())
            print(traceback.format_exc(),flush=True)
        finally:
            stop_children(child)
            time_exhausted()
            state.update(heartbeat=now(),finished_or_stopped_at=now())
            write(control/"status.json",state)
            print(f"{now()} {state['status']}: {state.get('reason','')}",flush=True)
        return 0 if state["status"]=="completed" else 1


def main(argv=None):
    for name in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","NUMEXPR_NUM_THREADS"):
        os.environ[name]="1"
    os.environ["MPLBACKEND"]="Agg"
    parser=argparse.ArgumentParser(description=__doc__);commands=parser.add_subparsers(dest="command",required=True)
    start=commands.add_parser("prepare",help="Freeze source/inputs and create a reviewable launch plan, without solvers")
    start.add_argument("--config",default=str(RUNTIME/"configs/research.toml"));start.add_argument("--workers",type=int,default=3)
    start.add_argument("--min-free-gib",type=float,default=15);start.add_argument("--stall-hours",type=float,default=12)
    start.add_argument("--profile",choices=("full","core"),default="full",help="Core runs construct/tables and A/B/C main comparisons only")
    start.add_argument("--max-run-hours",type=float,default=0,help="Cumulative active-time cap across resumes; 0 is unlimited, preserves saved work")
    start.add_argument("--output",required=True,type=Path)
    for command in ("start","resume","status","stop","check","_check","_supervise","_stage"):
        p=commands.add_parser(command);p.add_argument("output",type=Path,nargs="?")
        if command=="_stage":p.add_argument("--stage",required=True)
    args=parser.parse_args(argv)
    if args.command=="prepare":
        if not 0<args.min_free_gib<10000 or not 0<args.stall_hours<10000:parser.error("Resource limits must be positive finite values")
        if not math.isfinite(args.max_run_hours) or args.max_run_hours<0:parser.error("Runtime budget must be finite and nonnegative")
        prepare(args);return
    if args.output is None:
        latest=RUNTIME/"output/latest_background.json"
        if not latest.exists():parser.error("No latest run recorded; specify output directory")
        args.output=Path(read(latest)["output"])
    output=args.output.resolve()
    if args.command in ("start","resume"):launch(output)
    elif args.command=="status":print(json.dumps(read(output/"control/status.json"),ensure_ascii=False,indent=2))
    elif args.command=="stop":
        write(output/"control/stop.request",dict(requested_at=now()))
        print("Stop requested; supervisor will interrupt its own workers and keep results/checkpoints.")
    elif args.command=="check":
        spec=read(output/"control/launch.json")
        subprocess.run([spec["interpreter"],str(output/"snapshot/runtime/background_run.py"),"_check",str(output)],check=True)
    elif args.command=="_check":check_snapshot(output)
    elif args.command=="_supervise":return supervise(output)
    elif args.command=="_stage":stage_entry(output,args.stage)


if __name__=="__main__":raise SystemExit(main())
