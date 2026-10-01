from __future__ import annotations

"""Launch worker.py in private per-worker copies of the TaskA_v4 template.

python run_parallel.py --stage main            # 300 scenarios x 6 strategies, k=3
python run_parallel.py --stage sens --k 1      # 100 scenarios x OD_CEN/OD_IJSH at k=1
"""

import argparse
import os
import shutil
import subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
PYTHON = "/home/workenv/Power_Road_Simple/.venv/bin/python3"
TEMPLATE = "/home/workenv/TaskA_v4/template"
OD_ROOT = "/home/workenv/OD"
TABLES = "/home/workenv/TaskA_v4/results/tables.json"


def ensure_worker(worker_id: int) -> str:
    d = os.path.join(HERE, "workers", f"w{worker_id:02d}")
    if not os.path.exists(d):
        shutil.copytree(TEMPLATE, d, ignore=shutil.ignore_patterns("__pycache__", "results", "tapb_fallback.log", "full_log.txt"))
        os.makedirs(os.path.join(d, "results"), exist_ok=True)
    for f in ("task_c_od.py", "task_c_weighted_triangle.py"):
        shutil.copy(os.path.join(OD_ROOT, f), os.path.join(d, f))
    return d


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", choices=("main", "sens", "decomp", "shift"), default="main")
    ap.add_argument("--n-workers", type=int, default=12)
    ap.add_argument("--k", type=int, default=3)
    ap.add_argument("--n-scenarios", type=int, default=None)
    ap.add_argument("--variant", default=None)
    ap.add_argument("--tag", default="", help="result-folder suffix for an alternative critical set")
    ap.add_argument("--critical-buses", default="8,17")
    ap.add_argument("--shelter-nodes", default="24")
    ap.add_argument("--weight-buses", default=None)
    ap.add_argument("--pair-mode", default="depot_to_critical", choices=("depot_to_critical", "all_pairs"))
    ap.add_argument("--dir-offset", type=int, default=0, help="use worker dirs w{offset+i} (avoid dirs in use by another stage)")
    args = ap.parse_args()

    sfx = "" if args.pair_mode == "depot_to_critical" else "_allpairs"
    if args.tag:
        sfx = "_" + args.tag
    if args.stage == "main":
        n = args.n_scenarios or 300
        result_root = os.path.join(HERE, "results", "main_300" + sfx)
        strategies = "CEN,JSH,IJSH,OD_CEN,OD_JSH,OD_IJSH" if not sfx else ("OD_CEN,OD_IJSH" if not args.tag else "OD_CEN")
    elif args.stage == "shift":
        n = args.n_scenarios or 100
        result_root = os.path.join(HERE, "results", f"shift_{args.variant}" + sfx)
        strategies = "CEN,JSH,IJSH,OD_CEN,OD_IJSH" if not sfx else ("OD_CEN,OD_IJSH" if not args.tag else "OD_CEN")
    elif args.stage == "decomp":
        n = args.n_scenarios or 300
        result_root = os.path.join(HERE, "results", "decomp_300")
        strategies = "MIX_JSHroad_IJSHpower,MIX_IJSHroad_JSHpower"
    else:
        n = args.n_scenarios or 100
        result_root = os.path.join(HERE, "results", f"sens_k{args.k}" + sfx)
        strategies = "OD_CEN,OD_IJSH" if not sfx else "OD_CEN"
    os.makedirs(result_root, exist_ok=True)
    os.makedirs(os.path.join(HERE, "logs"), exist_ok=True)
    procs = []
    for w in range(args.n_workers):
        wd = ensure_worker(w + args.dir_offset)
        log = open(os.path.join(HERE, "logs", f"{args.stage}_{args.variant or ''}k{args.k}{sfx}_w{w:02d}.log"), "a", buffering=1, encoding="utf-8")
        cmd = [PYTHON, os.path.join(HERE, "worker.py"), "--worker-id", str(w), "--n-workers", str(args.n_workers),
               "--n-scenarios", str(n), "--result-root", os.path.join(result_root, f"worker_{w:02d}"),
               "--tables", TABLES, "--k", str(args.k), "--strategies", strategies, "--pair-mode", args.pair_mode, "--critical-buses", args.critical_buses, "--shelter-nodes", args.shelter_nodes] + (["--weight-buses", args.weight_buses] if args.weight_buses else []) + (["--variant", args.variant] if args.variant else [])
        p = subprocess.Popen(cmd, cwd=wd, stdout=log, stderr=subprocess.STDOUT,
                             env=dict(os.environ, PYTHONPATH=wd, PYTHONUNBUFFERED="1"))
        procs.append((w, p, log))
        print(f"launched {args.stage} k={args.k} worker {w:02d} pid={p.pid}", flush=True)
    bad = []
    for w, p, log in procs:
        rc = p.wait(); log.close()
        print(f"worker {w:02d} rc={rc}", flush=True)
        if rc:
            bad.append(w)
    if bad:
        raise SystemExit(f"failed workers: {bad}")
    print("stage complete", flush=True)


if __name__ == "__main__":
    main()
