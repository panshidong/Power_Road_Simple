#!/usr/bin/env python3
"""Explicit entry point. Importing this file does not start an experiment."""
import os

# Each worker gets one native BLAS thread; TAP-B has its own explicit budget.
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
             "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[name] = "1"
os.environ["MPLBACKEND"] = "Agg"

if __name__ == "__main__":
    from austin_runtime.cli import main
    main()
