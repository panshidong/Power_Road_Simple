# Task B Trade-off Workflow

## What This Adds
- Deterministic critical locations in `critical_location.json`
- CRI sensitivity over `(w_e, w_a)`
- Sensitivity over CRI threshold and critical-access threshold
- Critical-access restoration summaries based on `A_z(t)`
- Three decision-rule families:
  - single objective
  - weighted sum
  - guardrail
- One combined CSV and one efficiency-equity scatter plot

## Current Placeholder Assumptions
- Zone = node
- Bus-to-zone is derived from `bus_location.json`
- Crew dispatch now defaults to `new_bus_to_link.json` with `bus_dispatch_mode=link_only`
- Critical locations are fixed to nodes `[1, 24]`
- Accessibility recovery is `A_z(t) = min(1, TT0_z / TT_z(t))`
- A zone reaches critical access once `A_z(t)` reaches the configured threshold

## Main Entry
- Run `python3 runner.py`
- Or run `python3 tradeoff_runner.py`

## Outputs
- `tradeoff_summary.csv`
  - `row_kind=optimized` for the main optimized experiments
  - `row_kind=sensitivity` for re-evaluation under alternate weights/thresholds
  - `is_pareto=1` marks non-dominated solutions in the `triangle_area` vs `gini_restore` plane
- `tradeoff_scatter.png`
  - scatter of triangle area vs gini restore
  - dashed Pareto front
- `tradeoff_summary.md`
  - quick run summary and key settings

## Decision Rules Included
- Single-objective optimization:
  - `triangle`
  - `equity:var_restore`
  - `equity:gini_restore`
  - `equity:p90_restore`
  - `equity:maximin_time_avg_cri_loss` for maximin equity
  - `critical_access:p90_access_restore`
- Weighted-sum optimization:
  - `triangle + lambda * fairness_metric`
  - lambda values: `0.25, 0.5, 1.0`
  - includes the maximin-equity loss term `1 - min_time_avg_cri`
- Guardrail optimization:
  - minimize `triangle`
  - enforce fairness metric to stay below `0.95 * baseline`

## Where To Tweak Later
- `tradeoff_runner.py`
  - base sequence
  - SA iterations / cooling rate
  - dispatch mapping mode (`bus_dispatch_mode`)
  - dispatch source file (`bus_to_link_source`)
  - weight grids
  - threshold grids
  - guardrail factor
- `taskB_setup.py`
  - zone mapping
  - critical locations
