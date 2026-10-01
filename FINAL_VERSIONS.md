# Final versions (as of 2026-09-30)

This repository holds the code for three studies that share one coupled power–road recovery
simulator (IEEE 33-bus feeder, Sioux Falls road network, TAP-B traffic assignment). Each study's
final results were produced by a slightly different snapshot of that simulator, so each study has
its own branch. This file is identical on all four branches listed below.

Results, run logs, manuscripts, figures and the scripts that build or audit manuscripts are not in
the repository, because the repository is public. They stay in the author's local workspace
(`/home/workenv/`). The folder names below refer to that workspace.

## Which branch is final

| Study | Branch | Status | Local folder that produced the final results |
|---|---|---|---|
| Task A: pre-event criticality priorities (structural, service-based Shapley, access-augmented Shapley) | `TaskA_v4` | **final** | `TaskA_v4/`, plus supplementary checks in `TaskA_v5/` and `TaskA_v6/` |
| Task B: equity-aware recovery sequencing (CRI) | `TaskB_powerfix` | **final** | `TaskB_powerfix/` |
| Task C: origin–destination representation of road criticality | `TaskC_OD` | **final** | `OD/`, including `OD/taskc_offline/` |
| Task A, earlier design (scenario-specific S0/S1/S2 Shapley rankings) | `Criticality_use` | superseded by `TaskA_v4` | `Power_Road_Simple/` |

The other branches (`master`, `june26-submit`, `Equity`, `equity_modern_style`, `Criticality`,
`codex/fast100-specialized-batch`) predate the power-rule fix below and are kept for history.

## The power-rule fix (2026-09-03)

`power_util.get_functional_nodes` did not propagate a bus outage downstream on the feeder. It now
marks a bus unavailable when the bus or any upstream bus is broken, the same rule `delete_buses`
already used. `power_util.py` is byte-identical on `TaskA_v4`, `TaskB_powerfix` and `TaskC_OD`;
`Criticality_use` has the same logic with a different comment. Results produced before the fix are
"pre-fix" and must not be mixed with post-fix results in one table. The analysis of the bug is in
`backup_prefix_powerbug_20260903/BACKUP_RECORD.md` on `TaskC_OD`.

## Simulator differences between the final branches

| File | `TaskB_powerfix` | `TaskC_OD` (`task_c_*.py` runs) | `TaskA_v4` (also used by `OD/taskc_offline`) |
|---|---|---|---|
| `power_util.py` | fixed | fixed | fixed |
| `disaster.py` | original: road damage sampled over directed arcs | original | each physical road is one canonical undirected pair; both directions are damaged and repaired together |
| `run_tapb.py` | original: a TAP-B failure raises | original | if gap 1e-4 fails, retries with gaps 0.005, 0.01, 0.02, 0.05 and logs each fallback to `tapb_fallback.log` |
| `task_a_criticality.py` | original | original | structural road score uses the larger of the two directed betweenness values |
| `access_value.py` | absent | absent | new: depot-to-damaged-bus accessibility term of the access-augmented Shapley value |

Do not copy the `TaskA_v4` version of `disaster.py` into Task B. With the same seed it draws
different damage sets, and the Task B results would no longer reproduce.

## Layout of the final branches

- **`TaskA_v4`**: simulator at the repository root. `TaskA_v4/EXPERIMENT_SPEC.md` is the frozen
  experiment specification. Pipeline: `run_experiment.py` drives `construction_worker.py`
  (scenario Shapley scores), then `build_tables.py` (pre-event expected tables), then
  `evaluation_worker.py` and `robustness_worker.py`, then `analyze_results.py`.
  `TaskA_v5/` holds the representative-case and road-closure checks; `TaskA_v6/` holds the
  exploratory event-specific greedy comparator.
- **`TaskB_powerfix`**: simulator at the root (original plus the power-rule fix).
  `TaskB_powerfix/README.md` describes the parallel rerun of the main 100-scenario batch.
- **`TaskC_OD`**: Task C scripts (`task_c_*.py`) at the root. `taskc_offline/` re-evaluates the
  O–D road score against the `TaskA_v4` tables; it reads `TaskA_v4/template` and
  `TaskA_v4/results/tables.json` from the local workspace. `backup_prefix_powerbug_20260903/` is
  the pre-fix record.

## Scenario seeds

| Seed (seed0) | Used for |
|---|---|
| 20260325 | Task B main batch (pre-fix and post-fix runs use identical damage sets) |
| 20260402 | Task A earlier design (`Criticality_use`) |
| 20263001 | Task A v4 table-construction ensemble |
| 20264001 | Task A v4 evaluation ensemble (300 disruptions); also `OD/taskc_offline` |
| 20265001 / 20266001 / 20267001 / 20268001 | Task A v4 shifted ensembles: small / large / light / clustered |
| 20260701 | Task C main batch |
| 20260801 | Task C k-sensitivity |

## Rerunning

The task scripts were run inside the local workspace and contain absolute paths such as
`/home/workenv/TaskA_v4/template` and `/home/workenv/Power_Road_Simple/.venv/bin/python3`. The
simulator writes fixed scratch files (`s.txt`, `matrix0.bin`, `full_log.txt`,
`work/SiouxFalls_net1.txt`, `work/SiouxFalls_net2.txt`), so every parallel worker runs in a
private copy of the simulator, and two runs must never share a directory. To rerun from a fresh
clone, copy the branch root, including the `tap-b` submodule and its compiled binary, into
`<task folder>/template/` and adjust the absolute paths.
