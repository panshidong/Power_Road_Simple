# Power-functionality fix audit

Date: 2026-09-03

## Source reviewed

- External worktree: `/home/workenv/OD`
- Reviewed change: `OD/power_util.py`
- OD backup record: `OD/backup_prefix_powerbug_20260903/BACKUP_RECORD.md`
- OD recorded patch: `OD/backup_prefix_powerbug_20260903/fix.diff`

The OD worktree changes one pre-existing source file: `power_util.py`. The fix
replaces the traversal-order-dependent recursive implementation of
`get_functional_nodes` with the same directed downstream closure already used by
`delete_buses`.

## Merge decision

Only the power-functionality fix was merged. Task C/OD scripts and outputs were
not copied into the criticality worktree. The explanatory comment was shortened
to remove the obsolete “function is broken” note and to state the common outage
semantics directly.

## Regression checks

`task_a_checks.check_power_functionality_semantics` verifies:

- empty, root, intermediate, leaf, multiple-bus, and all-bus outages;
- every one of the 33 possible single-bus outages;
- equality between `get_functional_nodes(broken)` and the complement of
  `delete_buses(broken)`;
- the expected 24-bus downstream closure for a failure at bus 5.

The complete `task_a_checks.py` suite passed after the merge, including the
Shapley efficiency, strategy-generation, full-functionality, specialized-crew,
and dynamic TAP-B snapshot checks.

## Version preservation and rerun targets

- Pre-merge source and the v2 manuscript are preserved under
  `archive/powerfix_20260903_premerge/`.
- Original completed results remain in `results/task_a_criticality_final_100/`.
- The corrected rerun writes to `results/task_a_criticality_powerfix_100/`.
- The corrected stability check writes to
  `results/task_a_shapley_stability_powerfix_10/`.
- The corrected paper will be written as `Criticality_v3_powerfix.docx`.

The new `random_disasters.json` is byte-identical to the original 100-scenario
file, confirming that the rerun changes the power-performance rule while holding
the seeded damage sets and road capacity factors fixed.
