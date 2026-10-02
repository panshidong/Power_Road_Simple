"""Small, explicit progress messages for long destination-host operations."""
from __future__ import annotations

import contextlib
from datetime import datetime
import os
import threading
import time


def report(message):
    print(f"[{datetime.now().astimezone().isoformat(timespec='seconds')}] "
          f"pid={os.getpid()} {message}", flush=True)


@contextlib.contextmanager
def phase(label, *, enabled=True, interval=30):
    """Report entry, elapsed-time heartbeats, and exit; this is not a timeout."""
    if not enabled:
        yield
        return
    started = time.monotonic()
    stopped = threading.Event()
    report(f"START {label}")

    def heartbeat():
        while not stopped.wait(interval):
            report(f"RUNNING {label}; elapsed={time.monotonic()-started:.0f}s")

    thread = threading.Thread(target=heartbeat, daemon=True)
    thread.start()
    try:
        yield
    except BaseException:
        report(f"FAILED {label}; elapsed={time.monotonic()-started:.1f}s")
        raise
    else:
        report(f"DONE {label}; elapsed={time.monotonic()-started:.1f}s")
    finally:
        stopped.set()
        thread.join(timeout=1)
