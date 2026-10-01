# taskc_offline — O-D road criticality vs the pre-event Task A tables (OD paper v4)

Purpose: re-evaluate the Task C static O-D road score against the *final* Task A design
(pre-event tables CEN / JSH / IJSH from `/home/workenv/TaskA_v4/results/tables.json`,
checkpoint 400) on that study's own 300 evaluation disruptions (seeds 20264001+,
38 physical roads, corrected power rule, TAP-B gap-relaxation fallback).

Nothing in TaskA_v4, Power_Road_Simple or the OD worktree's own simulator files was
modified. Each worker directory `workers/wNN/` is a private copy of
`/home/workenv/TaskA_v4/template` plus `task_c_od.py` and `task_c_weighted_triangle.py`
copied from `/home/workenv/OD`.

Strategies per scenario (specialized crews → only within-trade order matters):
- `CEN`, `JSH`, `IJSH`: the Task A tables applied exactly as `TaskA_v4/evaluation_worker.py`
  (reproduce Codex's per-scenario losses with max |Δ| = 0).
- `OD_CEN`, `OD_JSH`, `OD_IJSH`: road block ordered by the static O-D score (k = 3, weights
  1.0/0.6/0.35, depot 1 → bus 8, bus 17 locations and shelter 24), power block ordered by the
  named table's power scores. Each `OD_X vs X` isolates the road-order change.
- `MIX_IJSHroad_JSHpower`, `MIX_JSHroad_IJSHpower` (stage `decomp`): attribute the IJSH-vs-JSH
  benefit to the road vs the power block.

Files
- `worker.py` — per-scenario simulation (resume-safe `rows_<scenario>.csv`), scenario manifest
  cross-checked against `TaskA_v4/results/evaluation/worker_00/evaluation_disasters.json`.
- `run_parallel.py --stage main|sens|decomp [--k K] [--dir-offset N]` — 12 private workers.
- `analyze.py` — merges rows → `results/merged/rows.csv`, `stats.json`, `SUMMARY.md`
  (paired bootstrap 20,000 resamples, trimmed 2.5 % per tail, medians, wins/losses/ties, exact sign test).
- `plot_figures.py` — `od4_*.png` in `/home/workenv/OD` for the paper builder.
- Results: `results/main_300/` (6 strategies × 300), `results/sens_k{1,2,5}/` (OD_CEN, OD_IJSH × first
  100 cases), `results/decomp_300/`, `results/merged/`.

Headline (300 cases, triangle loss): OD+CEN vs CEN −32.6 % (252/41/7); OD+JSH vs JSH −24.9 % (220/78/2);
OD+IJSH vs IJSH −6.5 % (172/114/14, sign p 0.0007); decomposition: IJSH's 19.9 % gain over JSH is
19.5 % road order, 1.2 % power order. k-sensitivity spread ≤ 3.8 %.

Paper: `/home/workenv/OD/build_od_paper.py` → `OD_Representation_Paper_Draft.docx` (v4); the previous
draft/builder are kept as `*_v3_onlineShapley.*`.

Shifted-damage ensembles (added later the same day): `run_parallel.py --stage shift --variant small|large|light|clustered`
→ `results/shift_<variant>/` (same definitions and seeds as `TaskA_v4/robustness_worker.py`; the three tables reproduce
Codex's robustness rows exactly), `analyze_shift.py` → `results/merged/shift_stats.json`, `SHIFT_SUMMARY.md`.
Headline: O-D vs access-augmented table −20.7 % (large), −14.0 % mean / 62-38 wins (clustered), −2.0 % (small), +3.0 % worse (light, no closures).

Alternative O-D configurations (v7): `--pair-mode all_pairs` (12 ordered trips among depot + facilities) →
`results/*_allpairs`; alternative critical sets via `--tag/--critical-buses/--shelter-nodes/--weight-buses`:
`roadcentral` (road nodes 6, 8, 16 = top betweenness), `central` (buses 2, 3 + node 6; power-centrality, not used in
the paper), `randA/B/C` (seeds 20260911–13), `size4/7/11` (seeds 20260914/17/21), `allnodes` (all 23 non-depot nodes).
`analyze.py --variant <tag>` merges each with the baselines → `results/merged/stats_<tag>.json`.
Finding: depot-origin trips carry the value (all-pairs ties the access-augmented table; larger or depot-adjacent sets do worse).
