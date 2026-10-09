Austin model-v2 completed core results — 2026-10-09

Run: core-v2-stallfix-20261008
Branch: austin-model-fixes
Numerical fingerprint: 21bd74ff9e413dc306d81dc4609e1574b63b18d6155732d54bb264735ef63dbe
Finished: 2026-10-09 03:42:39 America/Chicago
Machine runtime: 17 hours 12 minutes 30 seconds, with 12 workers and 2 TAP-B threads per call.

All 336 tasks completed: 72 construction, 72 A, 48 B, and 144 C.
All 264 evaluation trajectories completed every repair; unfinished repairs=0.
288 construction permutations and 504 simulated-annealing iterations completed.
0 failed tasks. Physical acceptance and final result audit passed.
The longest modeled recovery took 1627.013 minutes (about 27 hours 7 minutes).
This modeled recovery duration is separate from the machine runtime above.

DOWNLOAD ALL RESULTS

https://github.com/panshidong/Power_Road_Simple/raw/refs/heads/austin-model-fixes/reports/core_v2_results_20261009/Austin_model_v2_complete_results_20261009.tar.xz

Archive size: 12,904,124 bytes (12.31 MiB).
SHA256: 2c55b9bd4b1a22e8e4ff96ab458d44fd3186dbed3351fb10df31bb88f1204c45

This lossless archive contains 11,859 files. All 11,858 original
members were verified byte-for-byte against their source files. All 336 raw result
SHA256 hashes also match the completed-run acceptance record.

Extract into an empty directory in WSL/Linux:
  mkdir Austin_model_v2_results
  tar -xJf Austin_model_v2_complete_results_20261009.tar.xz -C Austin_model_v2_results
  cd Austin_model_v2_results
  sha256sum -c SHA256SUMS.txt

WHAT IS INCLUDED

original_run/results/       All 336 original JSON task results, including all
                            recovery events, repair dispatches and zone arrays.
original_run/analysis/      All original analysis tables, paired statistics and figures.
original_run/checkpoints/   All completed construction and optimization checkpoints.
original_run/plans/         Exact task plans.
original_run/control/       Supervisor logs, launch configuration, health review,
                            final acceptance and the actual completion auditor.
original_run/scratch/       All retained native request/result JSON and solver logs.
original_run/validation*    Physical AC/TAP-B validation and retained diagnostics.
code_snapshot/             Frozen Python runtime/configuration and exact native
                            executables, plus matching TAP-B source and regressions.
prepared_metadata/         Frozen catalog and load catalog.
solver_fix_reproduction/   Before/after TAP-B failure-case inputs and logs.
provenance/                Fix evidence, tests, completion record and input inventory.
SHA256SUMS.txt              Per-member integrity checks.

Browse the same analysis outputs directly on GitHub:
https://github.com/panshidong/Power_Road_Simple/tree/austin-model-fixes/reports/core_v2_results_20261009/analysis

Final acceptance and run configuration:
https://github.com/panshidong/Power_Road_Simple/tree/austin-model-fixes/reports/core_v2_results_20261009/run_metadata

The full public input/circuit copies and approximately 704 MiB of reusable physical
state caches remain preserved locally. They are not duplicated in this result
publication. provenance/input_inventory.json records the frozen input hashes;
code_snapshot/sources.lock.json identifies the public sources. Virtual environments,
lock files and Python bytecode are excluded. No result JSON, event series,
checkpoint, analysis output or retained solver diagnostic is omitted.

CODE AND INTERPRETATION

The archive's code_snapshot is the exact version associated with these results.
Repository source paths outside this result directory may not include all local
TAP-B patches; use the included matching source, executable hashes and fix evidence.
The publisher starts no solver and changes no original result or local cache.

The TAP-B adaptive bush-skipping starvation fix retains the original gap target
1e-4. The old failing Austin state now converges, with an independently checked
full-network gap 6.8456376851e-5 and node imbalance about 2e-6.

This is the core model-v2 matrix, not the full sensitivity/shift study. Model-v2 uses
healthy-base-case-referenced AC limits, local shedding, severity-based repair times,
damage-scaled crews and capped crew travel. There is no observation-window stop.
The automatically generated original SUMMARY.md has generic horizon-censoring text;
this run's actual completion audit reports all 264 trajectories uncensored and no
unfinished repairs. Disaster-size and other empirical calibration gaps remain.
These results should not be treated as calibrated real-world disaster predictions
or directly equated with the earlier model-v1 results or legacy article results.

archive_verification.json records the complete archive check. publish_manifest.json
records SHA256 and Git blob hashes for every published file except itself.
