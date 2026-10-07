from __future__ import annotations
import contextlib, csv, gzip, hashlib, importlib.metadata, json, math, os, sys, tempfile, tomllib
from pathlib import Path

RUNTIME = Path(__file__).resolve().parents[1]
AUSTIN = RUNTIME.parent


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def atomic_json(path, value):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(prefix=path.name+".", suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as f:
            json.dump(value, f, ensure_ascii=False, sort_keys=True, allow_nan=False)
            f.write("\n"); f.flush(); os.fsync(f.fileno())
        os.replace(name, path)
    finally:
        if os.path.exists(name): os.unlink(name)


def digest(obj):
    return hashlib.sha256(json.dumps(obj, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()).hexdigest()


def sha_file(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""): h.update(b)
    return h.hexdigest()


def rows(path):
    path = Path(path)
    with (gzip.open(path, "rt", encoding="utf-8", newline="") if path.suffix == ".gz" else path.open(encoding="utf-8", newline="")) as f:
        yield from csv.DictReader(f)


def write_rows(path, records):
    records = list(records); path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    with (gzip.open(path, "wt", encoding="utf-8", newline="") if path.suffix == ".gz" else path.open("w", encoding="utf-8", newline="")) as f:
        fields = list(dict.fromkeys(k for row in records for k in row))
        w = csv.DictWriter(f, fieldnames=fields); w.writeheader(); w.writerows(records)


def merge(a, b):
    out = dict(a)
    for k, v in b.items():
        out[k] = merge(a.get(k, {}), v) if isinstance(v, dict) else v
    return out


def config(path=None):
    with (RUNTIME/"configs/research.toml").open("rb") as f: cfg = tomllib.load(f)
    if path:
        with Path(path).open("rb") as f: cfg = merge(cfg, tomllib.load(f))
    if cfg["power"]["model"] != "tamu_regional_opendss": raise ValueError("Only the physical TAMU engine is implemented")
    recovery_attempts=cfg["power"].get("native_recovery_attempts",0)
    if type(recovery_attempts) is not int or recovery_attempts not in (0,1):raise ValueError("Native AC recovery attempts must be 0 or 1")
    if not 1 <= cfg["runtime"]["tapb_threads"] <= 64: raise ValueError("TAP-B threads must be in [1,64]")
    scales = cfg["power"]["load_scales"]
    if not scales or scales != sorted(set(scales), reverse=True) or not all(0 < x <= 1 for x in scales):
        raise ValueError("load_scales must be unique, descending, positive fractions")
    if cfg["runtime"]["workers"] < 0 or cfg["runtime"]["reserve_cpus"] < 0 or cfg["runtime"]["memory_gb_per_worker"] <= 0:
        raise ValueError("Invalid worker/resource budget")
    if cfg["runtime"]["jobs_per_worker_batch"] < 1: raise ValueError("jobs_per_worker_batch must be positive")
    if any(cfg["recovery"][k] != 1 for k in ("power_crews", "road_crews")):
        raise ValueError("The reproduced dispatch policy requires one crew per trade")
    if type(cfg["recovery"].get("stop_at_horizon",True)) is not bool:
        raise ValueError("recovery.stop_at_horizon must be a boolean")
    if any(type(cfg["recovery"][k]) not in (int,float) or not math.isfinite(cfg["recovery"][k]) or cfg["recovery"][k] <= 0
           for k in ("power_repair_minutes","road_repair_minutes","horizon_minutes","offroad_speed_kph")):
        raise ValueError("Durations and off-road speed must be finite and positive")
    if any(cfg["task_a"][k] < 1 for k in ("construction_scenarios","evaluation_scenarios","shift_scenarios","shapley_permutations")):
        raise ValueError("Task A sample counts must be positive")
    if cfg["task_b"]["scenarios"] < cfg["task_b"]["representative_scenario"] or cfg["task_b"]["representative_scenario"] < 1:
        raise ValueError("Representative Task B scenario must be inside its ensemble")
    if cfg["task_c"]["sensitivity_scenarios"] > cfg["task_a"]["evaluation_scenarios"]:
        raise ValueError("Task C sensitivity scenarios must be a subset of the main ensemble")
    if cfg["task_a"]["unseen_asset_policy"] not in ("error","centrality_fallback_recorded"):
        raise ValueError("Unknown unseen asset policy")
    if cfg["task_b"]["temperature"] <= 0 or not 0 < cfg["task_b"]["cooling"] <= 1:
        raise ValueError("Invalid annealing schedule")
    if min(cfg["equity"]["electric_weight"],cfg["equity"]["access_weight"]) < 0 or not abs(cfg["equity"]["electric_weight"]+cfg["equity"]["access_weight"]-1) < 1e-9:
        raise ValueError("CRI weights must be nonnegative and sum to one")
    if not 0 < cfg["power"]["voltage_min_pu"] < cfg["power"]["voltage_max_pu"]:
        raise ValueError("Invalid operating voltage interval")
    lo,hi=cfg["power"]["partial_rating_range"]
    if not 0 < lo <= hi <= 1 or not 0 <= cfg["power"]["partial_failure_probability"] <= 1:
        raise ValueError("Invalid partial power damage parameters")
    if cfg["task_b"]["sa_iterations"] < 0 or cfg["analysis"]["bootstrap_resamples"] < 1 or not 0 <= cfg["analysis"]["trim_each_tail"] < .5:
        raise ValueError("Invalid optimization/analysis count or trim fraction")
    if any(alpha < 0 for alpha in cfg["task_a"]["alpha_values"]) or not 1 <= cfg["task_a"]["representative_cases"] <= 3:
        raise ValueError("Invalid alpha or representative case count")
    if not cfg["task_c"]["path_weights"] or min(cfg["task_c"]["path_weights"]) < 0 or min([cfg["task_c"]["k"]]+cfg["task_c"]["k_sensitivity"]) < 1:
        raise ValueError("Invalid OD path weights/counts")
    if not 0 < cfg["traffic"]["signal_outage_capacity_factor"] <= 1 or cfg["traffic"]["time_to_minutes"] <= 0:
        raise ValueError("Invalid signal capacity factor or time conversion")
    return cfg


def local_path(cfg, name):
    path = Path(cfg["runtime"][name]).expanduser()
    return path.resolve() if path.is_absolute() else (RUNTIME/path).resolve()


def fingerprint(cfg, prepared):
    from .native_env import allocator_path
    allocator=allocator_path(cfg)
    # Worker count/output location do not alter experiments; physics and seeds do.
    numerical = {k:v for k,v in cfg.items() if k != "runtime"}
    code = {str(p.relative_to(RUNTIME)):sha_file(p) for p in sorted((RUNTIME/"austin_runtime").glob("*.py"))}
    versions={k:importlib.metadata.version(k) for k in ("numpy","scipy","networkx","OpenDSSDirect.py","dss-python","dss-python-backend","matplotlib")}
    return digest(dict(config=numerical,code=code,catalog=prepared,packages=versions,
        python=sys.version,native_tapb=sha_file(RUNTIME/"build/tap-b/bin/tap"),tapb_threads=cfg["runtime"]["tapb_threads"],
        native_allocator=sha_file(allocator) if allocator else None))


@contextlib.contextmanager
def file_lock(path, blocking=True):
    # Linux flock releases on process exit; unlike lock-directory schemes no stale ownership remains.
    import fcntl
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a+") as handle:
        flags = fcntl.LOCK_EX | (0 if blocking else fcntl.LOCK_NB)
        try: fcntl.flock(handle, flags)
        except BlockingIOError: raise RuntimeError(f"Another process owns {path}") from None
        try: yield
        finally: fcntl.flock(handle, fcntl.LOCK_UN)


def cached_json(path, compute):
    path = Path(path)
    if path.exists(): return read_json(path)
    with file_lock(path.with_suffix(".lock")):
        if path.exists(): return read_json(path)
        result = compute(); atomic_json(path, result); return result


def cpu_budget(cfg, requested=None):
    import math
    cpus = len(os.sched_getaffinity(0)) if hasattr(os,"sched_getaffinity") else (os.cpu_count() or 1)
    quota_file=Path("/sys/fs/cgroup/cpu.max")
    if quota_file.exists():
        quota,period=quota_file.read_text().split()[:2]
        if quota!="max":cpus=min(cpus,max(1,int(quota)//int(period)))
    budget = max(1, cpus-cfg["runtime"]["reserve_cpus"])
    threads = cfg["runtime"]["tapb_threads"]
    if threads > budget: raise ValueError("TAP-B thread count exceeds available CPU budget")
    available = None
    if Path("/proc/meminfo").exists():
        for line in Path("/proc/meminfo").read_text().splitlines():
            if line.startswith("MemAvailable:"): available = int(line.split()[1])*1024
    memory_limit=Path("/sys/fs/cgroup/memory.max");memory_current=Path("/sys/fs/cgroup/memory.current")
    if memory_limit.exists() and memory_current.exists():
        limit=memory_limit.read_text().strip()
        if limit!="max":
            remaining=max(0,int(limit)-int(memory_current.read_text()))
            available=min(available,remaining) if available is not None else remaining
    memory_workers = math.floor(available/(cfg["runtime"]["memory_gb_per_worker"]*2**30)) if available is not None else 1
    if memory_workers < 1:raise ValueError("Available RAM is below the configured per-worker estimate")
    automatic = min(max(1,budget//threads), memory_workers)
    if requested is not None and requested < 1: raise ValueError("--workers must be positive")
    workers = requested or cfg["runtime"]["workers"] or automatic
    if workers > automatic:
        raise ValueError(f"Requested {workers} workers, CPU/memory budget permits {automatic}; adjust explicit budgets if measured safe")
    regional_limit=min(budget,memory_workers,6)
    regional_requested=cfg["runtime"].get("validation_power_workers",0)
    if regional_requested<0 or regional_requested>regional_limit:
        raise ValueError(f"validation_power_workers exceeds CPU/RAM/region limit {regional_limit}")
    return dict(workers=workers,validation_power_workers=regional_requested or regional_limit,
                tapb_threads=threads,available_cpus=cpus,available_memory_bytes=available,
                memory_estimate_gb_per_worker=cfg["runtime"]["memory_gb_per_worker"])
