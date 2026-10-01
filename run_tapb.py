from __future__ import annotations

import os
import subprocess
from datetime import datetime, timezone
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

    # Algorithm B occasionally terminates on a bush-cycle degeneracy in very
    # heavily derated coalition states.  Preserve the production tolerance as
    # the first attempt, then use a deterministic, bounded relaxation only for
    # that failed state.  This is preferable to reusing a stale s.txt file or
    # silently dropping the disruption.  Every fallback is logged for QA.
    attempts = [str(gap)]
    if str(gap) == "0.0001":
        attempts.extend(["0.005", "0.01", "0.02", "0.05"])
    last_cp = None
    last_cmd = None
    for attempt_index, attempt_gap in enumerate(attempts):
        cmd = [tap_exe, attempt_gap, str(classes), netfile, tripfile]
        if attempt_index:
            cmd.append("1")
        cp = subprocess.run(
            cmd,
            cwd=cwd,
            stdout=(subprocess.DEVNULL if stdout_to_devnull else None),
            stderr=(subprocess.DEVNULL if stderr_to_devnull else None),
        )
        last_cp, last_cmd = cp, cmd
        if cp.returncode == 0:
            if attempt_index:
                log_root = cwd or os.getcwd()
                with open(os.path.join(log_root, "tapb_fallback.log"), "a", encoding="utf-8") as handle:
                    handle.write(
                        f"{datetime.now(timezone.utc).isoformat()}\t{netfile}\t"
                        f"primary_gap={gap}\tfallback_gap={attempt_gap}\n"
                    )
            return cp
    if check:
        raise RuntimeError(f"TAP-B failed after gaps {attempts!r} (returncode={last_cp.returncode}) cmd={last_cmd!r}")
    return last_cp
