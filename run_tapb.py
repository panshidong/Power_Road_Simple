from __future__ import annotations

import os
import subprocess
from typing import Optional


def run_tapb(
    netfile: str,
    tripfile: str,
    gap: str = "0.0001",
    classes: str = "1",
    *,
    tap_exe: str = "tap-b/bin/tap",
    cwd: Optional[str] = None,
    stdout_to_devnull: bool = True,
    stderr_to_devnull: bool = True,
    check: bool = True,
) -> subprocess.CompletedProcess:
    """
    Run TAP-B (C binary). Keeps IO-based workflow.

    Raises on:
      - missing executable / input files
      - nonzero return code (if check=True)
    """
    if not os.path.exists(tap_exe):
        raise FileNotFoundError(f"tap executable not found: {tap_exe!r}")
    if not os.path.exists(netfile):
        raise FileNotFoundError(f"netfile not found: {netfile!r}")
    if not os.path.exists(tripfile):
        raise FileNotFoundError(f"tripfile not found: {tripfile!r}")

    cmd = [tap_exe, str(gap), str(classes), netfile, tripfile]
    cp = subprocess.run(
        cmd,
        cwd=cwd,
        stdout=(subprocess.DEVNULL if stdout_to_devnull else None),
        stderr=(subprocess.DEVNULL if stderr_to_devnull else None),
    )
    if check and cp.returncode != 0:
        raise RuntimeError(f"TAP-B failed (returncode={cp.returncode}) cmd={cmd!r}")
    return cp
