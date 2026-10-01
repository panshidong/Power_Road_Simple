# TaskB_powerfix — parallel rerun of the equity (Task B) main batch with the corrected power rule

Created 2026-09-03. Self-contained; nothing here touches `Power_Road_Simple/`,
`/home/workenv/OD/`, or the original result folders (`/home/workenv/results/*`).

## Why

`power_util.get_functional_nodes` in the original code did not propagate a bus outage
downstream on the IEEE-33 feeder, so every resilience-triangle power term, the CRI
electric component, and every SA-optimized sequence in the equity paper were computed
under a rule different from the one the paper states. See
`/home/workenv/OD/backup_prefix_powerbug_20260903/BACKUP_RECORD.md` for the analysis.
This folder reruns the equity paper's main 100-scenario batch under the corrected rule.

## What is being rerun

Exactly the configuration behind the equity paper's Table 4 / `results/newruns`
(`runner_fast100.build_fast100_config()` on branch `Criticality_use` @ 9e09101):

- 100 random disruptions, `seed0 = 20260325`, 8–15 buses, 3–13 links, capacity drop 0.5–1.0
  (identical damage sets to the paper's batch — same generator, same seed).
- Experiments: single_triangle, single_gini_restore, single_maximin_time_avg_cri,
  single_p90_access_restore, weighted_gini_restore_l050,
  weighted_maximin_time_avg_cri_l050, guardrail_gini_restore, plus baseline_reference.
- SA: seed 0, max_iter 80, T0 1.4, alpha 0.97, swap neighbourhood.
- Specialized crews (1 power, 1 road), CRI weights (0.133, 0.867), thresholds 0.9,
  critical locations [1, 24], `bus_dispatch_mode = link_only`.

Code state: `template/` is a copy of the `Criticality_use` branch files (simulator,
Task B runners, data, `tap-b/` with compiled binary), with **one** change relative to
the branch: `power_util.get_functional_nodes` re-implemented as
"functional = all buses − (broken ∪ descendants of broken)" (same fix as in the
`/home/workenv/OD` worktree and, independently, in `Power_Road_Simple/`). No Task C
files are included. Not rerun here: the representative-case sensitivity
(`extra_sensitivity_runner.py`, Tables 5–7) — run it afterwards on the merged
scenario_084 run if needed.

## How it runs (and why it is parallel)

The simulator is single-threaded Python plus one TAP-B subprocess per repair event,
and it uses fixed relative scratch files (`work/SiouxFalls_net1/2.txt`, `s.txt`,
`bus_location.json`, `taskB_*.json`), so two runs in one directory corrupt each other.
`parallel_batch.py` therefore copies `template/` into `workers/wNN/`, starts one
`worker.py` per copy with `cwd` set to that copy, and assigns scenario i to worker
`(i-1) mod N`. Each worker writes `results/worker_NN/rows_<scenario_id>.csv` per finished
scenario (atomic rename), so re-running `parallel_batch.py` resumes. When all workers
exit, `merge_results.py` builds `results/merged/` with the same four artifacts the
original `batch_tradeoff_runner.py` produces (rows, manifest, aggregate, error-bar PNG,
summary).

```
/home/workenv/Power_Road_Simple/.venv/bin/python3 parallel_batch.py --n-workers 12
/home/workenv/Power_Road_Simple/.venv/bin/python3 merge_results.py   # re-merge any time
tail -f logs/w00.log                                                 # progress
```

## Layout

- `template/` — code + data snapshot (do not run anything in here)
- `workers/wNN/` — per-worker private copies (scratch; can be deleted after merge)
- `results/worker_NN/` — per-worker rows + per-scenario `scenario_runs/tradeoff_*/`
- `results/merged/` — final merged outputs
- `logs/wNN.log` — per-worker progress logs

## Comparing with the paper

Pre-fix reference numbers: `/home/workenv/results/newruns/aggregate_tradeoff_summary.csv`
(Table 4) and `/home/workenv/results/batch_tradeoff_20260326_151436/` (Section 5.1 text,
older simulator logic). Scenario ids and seeds are identical, so per-scenario paired
comparison against `newruns/scenario_tradeoff_rows.csv` is valid.
