"""One uncached regional AC state per OS process; invoked by PowerEngine only."""
from __future__ import annotations
import faulthandler
import os
from pathlib import Path
import resource
import sys
import traceback
from .common import atomic_json, digest, local_path, read_json
from .native_env import verify_loaded_allocator


def solve(request_path):
    request=read_json(request_path);folder=Path(request["directory"])
    with (folder/"native_crash.log").open("a") as log:
        log.write(f"pid={os.getpid()} region={request['region']} fatal_signals_only=True\n");log.flush()
        faulthandler.enable(file=log,all_threads=False)
        engine=None
        try:
            allocator=verify_loaded_allocator(request["config"])
            from .power import PowerEngine
            cfg=request["config"];catalog=read_json(local_path(cfg,"prepared")/"catalog.json")
            engine=PowerEngine(cfg,catalog,request["fingerprint"],folder/"dss",verbose=request["verbose"])
            # Direct call: never recursively create another isolated process.
            result=engine._operate(request["region"],request["damage"])
            engine.close();engine=None
            atomic_json(folder/"result.json",dict(request_digest=digest(request),fingerprint=request["fingerprint"],
                        pid=os.getpid(),result=result,native_allocator_sha256=allocator,
                        peak_rss_mib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024))
        except BaseException:
            atomic_json(folder/"python_error.json",dict(error=traceback.format_exc()))
            raise
        finally:
            if engine is not None:engine.close(cancel=True)
            faulthandler.disable()


if __name__=="__main__":solve(Path(sys.argv[1]))
