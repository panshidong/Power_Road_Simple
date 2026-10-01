# Criticality paper v4 — frozen experiment specification

## Scope

The experiment evaluates reusable, pre-event criticality rankings for the
coupled power-road recovery problem.  It does not compare subsystem models or
attempt to identify which dependency channel is more important.  Every
strategy is executed by the same fully coupled event-driven simulator with one
specialized power crew and one specialized road crew.

## Main strategies

- **CEN**: structural criticality (radial downstream power impact and weighted
  road-edge betweenness).
- **JSH**: an offline expected joint Shapley table based on the coupled service
  value `0.5 * F_power + 0.5 * F_road`.
- **IJSH**: an offline expected interdependency-aware joint Shapley table based
  on `(0.5 * F_power + 0.5 * F_road + alpha * A) / (1 + alpha)`, where `A`
  measures normalized depot accessibility to the initially disrupted power
  facilities under the coalition state.  `alpha = 1` is the prespecified
  equal-range setting; other values are sensitivity checks, not optimization.

The global joint scores are consumed as trade-compatible priority orders by
the specialized crews.  Thus cross-system marginal contributions can change
the within-power and within-road priorities, while trade eligibility remains
realistic.

## Hypotheses

- **H1:** JSH produces lower system resilience-loss area than CEN.
- **H2:** IJSH produces lower system resilience-loss area than JSH.

## Monte Carlo design

- Table-construction ensemble: independent disruptions drawn from the stated
  hazard model; used only to estimate `E[phi_i(D) | i is damaged]`.
- Independent evaluation ensemble: new draws from the same model; used only to
  evaluate the frozen tables.
- Shifted-distribution ensembles: secondary robustness evidence.
- The samples are Monte Carlo integration/evaluation ensembles, not machine
  learning training and test data.

## Road semantics

A road repair is one canonical physical pair `(min(u,v), max(u,v))`.  Damage,
capacity adjustment, coalition repair, crew dispatch, and completion all refer
to that same physical component; both directed travel arcs are affected.

## Claims

All reported performance claims are conditional on the modeled networks,
hazard distributions, dependencies, and recovery resources.  The paper argues
that the strategy is promising and adaptable, not universally dominant or
optimal.
