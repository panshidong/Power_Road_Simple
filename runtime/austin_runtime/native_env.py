"""Pinned, optional allocator for isolated AC children only."""
from pathlib import Path
import os
from .common import RUNTIME, sha_file


def allocator_path(cfg):
    power=cfg.get("power",{})
    name=power.get("native_allocator_path","")
    expected=power.get("native_allocator_sha256","")
    if not name and not expected:return None
    if not power.get("isolate_states",False):
        raise ValueError("A native allocator requires isolated AC states")
    if (not isinstance(name,str) or not name or any(c.isspace() or c==":" for c in name)
        or not isinstance(expected,str) or len(expected)!=64
        or any(c not in "0123456789abcdef" for c in expected)):
        raise ValueError("Native allocator requires a loader-safe path and pinned SHA256")
    path=Path(name)
    if not path.is_absolute():path=RUNTIME/path
    path=path.resolve(strict=True)
    if sha_file(path)!=expected:raise ValueError("Native allocator SHA256 mismatch")
    return path


def child_environment(cfg, env):
    env=dict(env)
    path=allocator_path(cfg)
    if path:
        # Do not silently combine an untracked preload with the pinned library.
        inherited=env.get("LD_PRELOAD","")
        if inherited and inherited!=str(path):raise ValueError("Conflicting inherited LD_PRELOAD")
        env["LD_PRELOAD"]=str(path)
        env["MALLOC_CONF"]="background_thread:false"
    return env


def verify_loaded_allocator(cfg):
    path=allocator_path(cfg)
    if path is None:return None
    # The dynamic loader can warn and ignore a bad preload yet exit zero. Require
    # the exact checked object in this process's mappings before solving.
    mapped=[]
    for line in Path("/proc/self/maps").read_text().splitlines():
        fields=line.split(maxsplit=5)
        if len(fields)==6 and fields[5].startswith("/"):mapped.append(fields[5])
    if str(path) not in mapped or os.environ.get("LD_PRELOAD")!=str(path):
        raise RuntimeError("Pinned native allocator was not loaded; refusing AC solve")
    return cfg["power"]["native_allocator_sha256"]
