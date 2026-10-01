from __future__ import annotations

"""
Driver: launch N isolated workers (each a private copy of template/) that split
the 100 fast100 scenarios, wait for all, then merge. Re-running is safe: workers
skip finished scenarios.

  python parallel_batch.py --n-workers 12
"""

import argparse
import os
import shutil
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
PY = "/home/workenv/Power_Road_Simple/.venv/bin/python3"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-workers", type=int, default=12)
    args = ap.parse_args()

    template = os.path.join(HERE, "template")
    workers_root = os.path.join(HERE, "workers")
    results_root = os.path.join(HERE, "results")
    logs_root = os.path.join(HERE, "logs")
    for d in (workers_root, results_root, logs_root):
        os.makedirs(d, exist_ok=True)

    procs = []
    for k in range(args.n_workers):
        wdir = os.path.join(workers_root, f"w{k:02d}")
        if not os.path.exists(wdir):
            shutil.copytree(template, wdir, symlinks=False)
        log = open(os.path.join(logs_root, f"w{k:02d}.log"), "a", buffering=1)
        cmd = [PY, os.path.join(HERE, "worker.py"), "--worker-id", str(k), "--n-workers", str(args.n_workers),
               "--result-root", os.path.join(results_root, f"worker_{k:02d}")]
        env = dict(os.environ, PYTHONPATH=wdir, PYTHONUNBUFFERED="1")
        p = subprocess.Popen(cmd, cwd=wdir, stdout=log, stderr=subprocess.STDOUT, env=env)
        procs.append((k, p, log))
        print(f"launched worker {k} pid {p.pid} in {wdir}", flush=True)

    t0 = time.time()
    failed = []
    for k, p, log in procs:
        rc = p.wait()
        log.close()
        print(f"worker {k} exited rc={rc} after {(time.time()-t0)/60:.1f} min", flush=True)
        if rc != 0:
            failed.append(k)
    if failed:
        print(f"FAILED workers: {failed} — check logs/ and re-run parallel_batch.py to resume", flush=True)
        sys.exit(1)

    subprocess.run([PY, os.path.join(HERE, "merge_results.py")], check=True, cwd=HERE)
    print("ALL WORKERS DONE AND MERGED", flush=True)


if __name__ == "__main__":
    main()
