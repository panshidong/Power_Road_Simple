# AI Coding Agent Instructions for Power-Road Resilience Analysis

## Project Overview
This codebase analyzes resilience of interdependent power (IEEE 33-bus) and road (SiouxFalls) networks under disruptions. It optimizes restoration sequences using simulated annealing, considering crew scheduling, interdependencies, and equity metrics.

## Key Components
- **resilience_measurement.py**: Core evaluation engine integrating power/road simulations
- **simulated_annealing.py**: Generic SA optimizer for permutation sequences
- **scheduler.py**: Crew-based restoration scheduling with travel times from TAP-B
- **interdependency.py**: Maps power failures to road capacity reductions
- **crews.py**: Manages specialized/multifunction repair crews
- **tap-b/**: C-based traffic assignment solver (build with `cd tap-b && make`)

## Data Flow
1. Disruption: Broken buses/links → power_util.delete_buses() → unfunctional nodes
2. Interdependency: unfunctional_nodes → interdependency.power_to_road() → road_util.capacity_adjustment()
3. Traffic: Modified network → run_tapb() → s.txt (flows/travel times)
4. Resilience: eval_power_resilience() + eval_road_resilience() → functionality metrics
5. Scheduling: sequence + crews → timeline with travel/service times

## Common Patterns
- **Assets**: `int` for power buses (1-33), `(int, int)` for road links (e.g., `(1, 3)`)
- **Sequences**: Lists of assets, e.g., `[11, (1, 3), 15]` - restoration order
- **Objectives**: `"triangle"` (resilience loss area), `"equity:gini"` (inequality metric)
- **Crew Modes**: `"specialized"` (separate power/road crews), `"multifunction"` (combined skills)
- **Mappings**: `bus_to_link.json` (bus→link), `bus_location.json` (bus→node)

## Workflows
- **Build TAP-B**: `cd tap-b && make` (required for traffic simulation)
- **Run SA Optimization**: Use `optimize_sequence_sa()` from resilience_measurement.py
- **Evaluate Single Sequence**: Call `evaluate_with_crews()` from scheduler.py
- **Debug**: Check `s.txt` for TAP-B output, results/*/ for logs/reports

## Examples
- Basic SA run: `SAConfig(seed=0, max_iter=30, T0=1.0, alpha=0.95, neighbor="swap")`
- Crew pool: `make_crews(mode="multifunction", multifunction_crews=1)`
- Objective extraction: `equity["gini"]` for Gini coefficient from equity metrics

## Conventions
- Use `strict=True` in evaluations to catch missing files/mappings
- Results saved to `results/` with timestamps (e.g., `opt_sa_20260304_104633/`)
- Travel times computed from TAP-B's `s.txt` using shortest paths
- Functionality: 0.0 (fully disrupted) to 1.0 (fully functional)</content>
<parameter name="filePath">/home/workenv/Power_Road_Simple/.github/copilot-instructions.md