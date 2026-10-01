from __future__ import annotations

"""Parallel driver for the v4 table-construction and evaluation ensembles."""

import argparse
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PYTHON = "/home/workenv/Power_Road_Simple/.venv/bin/python3"


def ensure_workers(n_workers: int) -> None:
    os.makedirs(os.path.join(HERE, "workers"), exist_ok=True)
    for worker_id in range(n_workers):
        worker_dir = os.path.join(HERE, "workers", f"w{worker_id:02d}")
        if not os.path.exists(worker_dir):
            shutil.copytree(os.path.join(HERE, "template"), worker_dir)


def run_workers(script: str, stage: str, n_workers: int, n_scenarios: int, extra=()) -> None:
    ensure_workers(n_workers)
    result_root = os.path.join(HERE, "results", stage)
    log_root = os.path.join(HERE, "logs")
    os.makedirs(result_root, exist_ok=True)
    os.makedirs(log_root, exist_ok=True)
    jobs = []
    for worker_id in range(n_workers):
        worker_dir = os.path.join(HERE, "workers", f"w{worker_id:02d}")
        worker_results = os.path.join(result_root, f"worker_{worker_id:02d}")
        log_path = os.path.join(log_root, f"{stage}_w{worker_id:02d}.log")
        log = open(log_path, "a", buffering=1, encoding="utf-8")
        command = [
            PYTHON,
            os.path.join(HERE, script),
            "--worker-id", str(worker_id),
            "--n-workers", str(n_workers),
            "--n-scenarios", str(n_scenarios),
            "--result-root", worker_results,
            *extra,
        ]
        process = subprocess.Popen(
            command,
            cwd=worker_dir,
            stdout=log,
            stderr=subprocess.STDOUT,
            env=dict(os.environ, PYTHONPATH=worker_dir, PYTHONUNBUFFERED="1"),
        )
        jobs.append((worker_id, process, log))
        print(f"launched {stage} worker {worker_id:02d}, pid={process.pid}", flush=True)

    failures = []
    for worker_id, process, log in jobs:
        return_code = process.wait()
        log.close()
        print(f"{stage} worker {worker_id:02d} rc={return_code}", flush=True)
        if return_code:
            failures.append(worker_id)
    if failures:
        raise SystemExit(f"{stage} failed workers: {failures}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", choices=("construct", "tables", "evaluate", "robustness", "all"), default="all")
    ap.add_argument("--n-workers", type=int, default=12)
    ap.add_argument("--n-construction", type=int, default=400)
    ap.add_argument("--n-evaluation", type=int, default=300)
    ap.add_argument("--n-robustness", type=int, default=100, help="scenarios per shifted distribution")
    args = ap.parse_args()

    table_path = os.path.join(HERE, "results", "tables.json")
    if args.stage in ("construct", "all"):
        run_workers("construction_worker.py", "construction", args.n_workers, args.n_construction)
    if args.stage in ("tables", "all"):
        subprocess.run(
            [
                PYTHON,
                os.path.join(HERE, "build_tables.py"),
                "--construction-root", os.path.join(HERE, "results", "construction"),
                "--output", table_path,
                "--expected", str(args.n_construction),
            ],
            check=True,
            cwd=HERE,
        )
    if args.stage in ("evaluate", "all"):
        run_workers(
            "evaluation_worker.py",
            "evaluation",
            args.n_workers,
            args.n_evaluation,
            extra=("--tables", table_path, "--checkpoint", str(args.n_construction)),
        )
    if args.stage in ("robustness", "all"):
        run_workers(
            "robustness_worker.py",
            "robustness",
            args.n_workers,
            args.n_robustness,
            extra=("--tables", table_path, "--checkpoint", str(args.n_construction)),
        )
    print("v4 requested stage(s) complete", flush=True)


if __name__ == "__main__":
    main()
