# Backup record — pre-fix Task C state (power-performance rule bug)

Created: 2026-09-03, on git worktree `/home/workenv/OD`, branch `TaskC_OD`
(HEAD 9e09101 "Add Task A criticality workflow"; all Task C files untracked — see
`git_state.txt`). Nothing in `/home/workenv/Power_Road_Simple` (Codex's working
directory, branch `Criticality_use`) was touched.

## Why this backup exists

`power_util.get_functional_nodes()` — the function `resilience_measurement.eval_power_resilience`
uses to compute power performance P_power(t) = |functional buses| / 33 — did not implement the
model stated in the papers ("a bus is functional iff it and all upstream buses are functional").
It recursed over *children* (parent→children dictionary) with a `visited` set shared across
calls, so an outage was **not propagated downstream**: only the damaged buses, the feeder root
(bus 1) and a small, traversal-order-dependent set of ancestors were counted non-functional.
The file's own comment said "this function is broken, use the below one"; `delete_buses()`
(downstream propagation) existed but `eval_power_resilience` still called the broken one.

Measured on the 100 main-study initial damage sets (pre-fix results):

| quantity | old rule | downstream-propagation rule |
|---|---|---|
| mean initial P_power over 100 scenarios | 0.453 | 0.090 |
| broken {5}: non-functional buses | [1, 5] | 5 and 23 descendants |
| broken {2, 19} | [1, 2, 19] | 32 buses |
| broken {3, 30} | [1, 3, 4, 30] | 29 buses |

Scope of the old rule (all identical in `master`, `june26-submit`, `equity_modern_style`,
`Criticality_use`; unchanged since First commit 2024-08-25):
- `resilience_measurement.eval_power_resilience` → power term of every resilience triangle
  (Task A, B, C; SA objective).
- `resilience_measurement.run_model_multi` → `compute_E_by_zone` → CRI electric component.
- `task_a_criticality.FullFunctionalityValue` → S1/S2 Shapley value functions.
- `task_c_weighted_triangle.weighted_power_func` (Task C supplemental metric).
Not affected: road-side propagation of power outages (`interdependency.power_to_road` uses
`delete_buses`, which does propagate), P_road, crew dispatch, and the S0/S3 power-side
ranking score (`power_radial_service_centrality`, which uses `delete_buses`).

## The fix applied after this backup (only in the `/home/workenv/OD` worktree)

`power_util.get_functional_nodes(broken_nodes)` re-implemented as: unfunctional = broken ∪ all
descendants of broken (iterative traversal over `connections`); functional = all 33 − unfunctional.
Signature and return type unchanged. `delete_buses()` unchanged (it ignores the value it gets
from `get_functional_nodes`). The original file is preserved verbatim in
`original_code/power_util.py`; the applied diff is recorded in `fix.diff` (added after the edit).

This is the only change to any pre-existing file on the branch. Task A/B results in the main
repository were NOT rerun and still reflect the old rule.

## Contents

- `results/task_c_main_100/` — pre-fix main study: 100 scenarios (seed0 20260701, 8–15 buses,
  3–13 links, capacity drop 0.5–1.0), S0–S3, standard + weighted triangle + CRI metrics.
  Files: `scenario_manifest.csv`, `scenario_strategy_rows.csv` (400 rows),
  `criticality_rankings.csv`, `aggregate_task_c_summary.csv`, `task_c_summary.md`,
  `random_disasters.json`, `scenario_runs/` (baseline_s.txt and state snapshots per run).
- `results/task_c_sensitivity/` — pre-fix k-sweep (k=1,2,3,5; 15 scenarios, seed0 20260801;
  S3 only).
- `analysis_prefix.txt` — output of `task_c_analysis.py` on the pre-fix main study (paired stats).
- `paper/OD_Representation_Paper_Draft.docx` — v2 draft built from pre-fix results (the one sent
  to the user on 2026-09-03), plus the v1 draft (`_OLD`) and the three figures.
- `code/` — snapshot of all Task C scripts and `TASKC_ASSUMPTIONS.md` as of the backup.
- `original_code/` — verbatim `power_util.py`, `resilience_measurement.py`, `task_a_criticality.py`
  before the fix.
- `md5sums.txt`, `git_state.txt`.

## Key pre-fix numbers (for the "what changes" comparison)

Main study, 100 scenarios, mean ± 95% CI:

| strategy | standard triangle | weighted triangle | Gini | min time-avg CRI | P90 access |
|---|---|---|---|---|---|
| S0 | 2033.6 ± 720.2 | 3216.6 ± 1185.1 | 0.185 | 0.291 | 1685.2 |
| S1 | 1524.7 ± 199.5 | 2430.9 ± 325.5 | 0.168 | 0.260 | 1382.4 |
| S2 | 1558.5 ± 457.9 | 2413.1 ± 741.1 | 0.187 | 0.329 | 1377.8 |
| S3 | 1268.9 ± 373.1 | 1902.3 ± 583.7 | 0.220 | 0.354 | 1187.0 |

Paired S3 − baseline, standard triangle: vs S0 −764.7 ± 412.9 (−37.6 %, sig, 89/9/2);
vs S1 −255.8 ± 332.7 (−16.8 %, n.s., 89/11/0); vs S2 −289.6 ± 148.1 (−18.6 %, sig, 70/30/0).
Weighted triangle: vs S0 −1314.3 ± 699.2 (sig); vs S1 −528.5 ± 519.5 (sig); vs S2 −510.8 ± 248.4 (sig).
Static-score coverage: 290/793 damaged links nonzero (36.6 %); S3 road order == S0 in 2/100.
k-sweep (15 scen.): triangle 1695.7 / 1712.4 / 1711.3 / 1676.8 for k = 1/2/3/5.
