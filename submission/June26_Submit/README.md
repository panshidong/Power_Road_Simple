# June26 Submit

This folder collects the code-linked artifacts for the June 26 submission package.

## Figure And Table Sources

### 100-scenario aggregate trade-off figure

- Repository copy: `figures/aggregate_tradeoff_errorbars.png`
- Full source data copy: `figures/aggregate_tradeoff_summary.csv`
- Plotted data copy, excluding the variance objective: `figures/aggregate_tradeoff_summary_plotted.csv`
- Original run folder: `/home/workenv/results/batch_tradeoff_20260326_151436`
- Original source CSV: `/home/workenv/results/batch_tradeoff_20260326_151436/aggregate_tradeoff_summary.csv`
- Original figure path: `/home/workenv/results/batch_tradeoff_20260326_151436/aggregate_tradeoff_errorbars.png`
- Scenario count: 100 random disaster scenarios.
- Figure content: aggregate mean triangle area versus mean Gini restoration, with 95% confidence intervals.
- Plotting note: the figure was regenerated at 600 DPI with larger fonts and descriptive labels. The variance objective point was intentionally removed from the displayed figure; use `aggregate_tradeoff_summary_plotted.csv` to reproduce the plotted rows exactly.

### Current-code representative special-case trade-off figure

- Repository copy: `figures/current_code_special_case_scatter.png`
- Source data copy: `special_case/current_code_tradeoff_summary.csv`
- Source summary copy: `special_case/current_code_tradeoff_summary.md`
- Source notes copy: `special_case/current_code_notes.md`
- Original run folder: `/home/workenv/results/batch_tradeoff_20260326_151436/scenario_runs/tradeoff_20260326_180346/current_code_special_case_20260604/tradeoff_20260604_141137`
- Original figure path: `/home/workenv/results/batch_tradeoff_20260326_151436/scenario_runs/tradeoff_20260326_180346/current_code_special_case_20260604/tradeoff_20260604_141137/tradeoff_scatter.png`
- Scenario: `scenario_084`, seed `20260408`.
- Disruption: 9 failed power buses and 10 damaged road links.
- CRI weights: `(w_e, w_a) = (0.133, 0.867)`.
- Weighted lambda in the main weighted runs: `0.5`.
- Figure content: current-code special-case resilience-equity scatter plot, triangle area versus Gini restoration.
- Plotting note: this is the preferred special-case figure because it was regenerated with the current dynamic evaluator. The variance objective point was intentionally removed from the displayed figure.

### Sensitivity analysis document

- Repository copy: `sensitivity/Sensitivity_updated.md`
- Word document copy: `sensitivity/Sensitivity_updated.docx`
- Source extra-sensitivity summary: `sensitivity/extra_sensitivity_summary.md`
- Source weighted lambda rows: `sensitivity/weighted_lambda_rows.csv`
- Original extra-sensitivity run folder: `/home/workenv/results/batch_tradeoff_20260326_151436/scenario_runs/tradeoff_20260326_180346/extra_sensitivity_20260604_123941`
- Original document paths:
  - `/home/workenv/results/Sensitivity_updated_20260603.md`
  - `/home/workenv/results/Sensitivity_updated_20260603.docx`
- Sensitivity settings:
  - CRI weight sensitivity: `(0.3, 0.7)` and `(0.7, 0.3)`.
  - Threshold sensitivity: `0.85` and `0.95`.
  - Weighted lambda sensitivity around the original `lambda = 0.5`: `0.25` and `0.75`.

## Important Consistency Note

The original special-case scatter plot under `/home/workenv/results/batch_tradeoff_20260326_151436/scenario_runs/tradeoff_20260326_180346/tradeoff_scatter.png` should not be mixed with the updated sensitivity values. The original run artifacts contain a single `dispatch_s.txt`, while the current-code rerun contains per-state `state_snapshots/`, reflecting the later dynamic simulation correction in which TAP-B travel times are updated after repair events.

For the representative special case, use the current-code rerun listed above.
