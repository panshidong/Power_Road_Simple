"""Explicit commands only: importing this module starts no work."""
from __future__ import annotations
import argparse, json, os
from .common import config, file_lock, fingerprint, local_path, read_json
from .plan import STAGES, jobs


def main(argv=None):
    for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        os.environ[name] = "1"
    parser = argparse.ArgumentParser(description="Austin physical AC power / TAP-B recovery experiments")
    parser.add_argument("--config", help="TOML overlay; default is runtime/configs/research.toml")
    parser.add_argument("--output", help="Absolute path or path relative to runtime/")
    commands = parser.add_subparsers(dest="command", required=True)
    for name in ("plan", "doctor", "prepare", "validate", "analyze"):
        commands.add_parser(name)
    run_parser = commands.add_parser("run")
    run_parser.add_argument("--stage", choices=["all"]+STAGES, default="all")
    run_parser.add_argument("--workers", type=int)
    args = parser.parse_args(argv)
    cfg = config(args.config)
    if args.output: cfg["runtime"]["output"] = args.output
    if args.command == "plan":
        counts = {stage:len(jobs(cfg,stage)) for stage in STAGES if stage not in ("tables","analyze","a-cases")}
        counts["a-cases"] = cfg["task_a"]["representative_cases"]*len(cfg["task_a"]["closure_penalties"])*3
        result = dict(stages=counts,total_jobs=sum(counts.values()),output=str(local_path(cfg,"output")),
            construction_permutations=cfg["task_a"]["construction_scenarios"]*cfg["task_a"]["shapley_permutations"],
            note="Counts only. Does not prepare data, import native engines, or run experiments. Cache misses determine native solve count.")
    elif args.command == "doctor":
        from .validation import doctor
        result = doctor(cfg)
    elif args.command == "prepare":
        from .prepare import prepare
        result = prepare(cfg)
    elif args.command == "validate":
        from .validation import validate
        output=local_path(cfg,"output")
        with file_lock(output/"run.lock",blocking=False):
            manifest=output/"run_manifest.json"
            if manifest.exists() and read_json(manifest)["fingerprint"]!=fingerprint(cfg,read_json(local_path(cfg,"prepared")/"catalog.json")):
                raise ValueError("Existing experiment has another fingerprint; use a new output")
            result = validate(cfg)
    else:
        from .workflow import run
        result = run(cfg, "analyze" if args.command == "analyze" else args.stage,
                     None if args.command == "analyze" else args.workers)
    print(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
