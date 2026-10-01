# Post-fix change list (pre-fix → post-fix)

Pre-fix scenarios: 100; post-fix scenarios: 100 (same seed0 → identical damage sets).

## 0. Did the repair sequences themselves change?

| strategy | scenarios where (power seq, road seq) differs old vs new |
|---|---|
| S0 | 0 / 100 |
| S1 | 100 / 100 |
| S2 | 100 / 100 |
| S3 | 0 / 100 |

Expected: S0 and S3 unchanged (their scores do not use the fixed function); S1/S2 may change (their Shapley value functions use power performance).

## 1. Marginal means (mean ± 95% CI)

### triangle_area

| strategy | pre-fix | post-fix |
|---|---|---|
| S0 | 2033.6 ± 720.2 | 2376.5 ± 831.8 |
| S1 | 1524.7 ± 199.5 | 1703.6 ± 472.5 |
| S2 | 1558.5 ± 457.9 | 1779.9 ± 527.9 |
| S3 | 1268.9 ± 373.1 | 1539.3 ± 453.3 |

Best-to-worst order: pre-fix S3 < S1 < S2 < S0 → post-fix S3 < S1 < S2 < S0

### weighted_triangle_area

| strategy | pre-fix | post-fix |
|---|---|---|
| S0 | 3216.6 ± 1185.1 | 3576.2 ± 1294.7 |
| S1 | 2430.9 ± 325.5 | 2476.3 ± 697.3 |
| S2 | 2413.1 ± 741.1 | 2596.1 ± 800.5 |
| S3 | 1902.3 ± 583.7 | 2183.8 ± 658.1 |

Best-to-worst order: pre-fix S3 < S2 < S1 < S0 → post-fix S3 < S1 < S2 < S0

### gini_restore

| strategy | pre-fix | post-fix |
|---|---|---|
| S0 | 0.185 ± 0.012 | 0.160 ± 0.008 |
| S1 | 0.168 ± 0.010 | 0.160 ± 0.008 |
| S2 | 0.187 ± 0.012 | 0.159 ± 0.009 |
| S3 | 0.220 ± 0.010 | 0.176 ± 0.008 |

Best-to-worst order: pre-fix S1 < S0 < S2 < S3 → post-fix S2 < S0 < S1 < S3

### min_time_avg_cri

| strategy | pre-fix | post-fix |
|---|---|---|
| S0 | 0.291 ± 0.025 | 0.289 ± 0.025 |
| S1 | 0.260 ± 0.017 | 0.336 ± 0.021 |
| S2 | 0.329 ± 0.021 | 0.337 ± 0.022 |
| S3 | 0.354 ± 0.020 | 0.358 ± 0.019 |

Best-to-worst order: pre-fix S3 < S2 < S0 < S1 → post-fix S3 < S2 < S1 < S0

### p90_access_restore

| strategy | pre-fix | post-fix |
|---|---|---|
| S0 | 1685.2 ± 512.9 | 1685.2 ± 512.9 |
| S1 | 1382.4 ± 147.9 | 1293.0 ± 304.8 |
| S2 | 1377.8 ± 317.9 | 1321.0 ± 309.5 |
| S3 | 1187.0 ± 302.1 | 1187.0 ± 302.1 |

Best-to-worst order: pre-fix S3 < S2 < S1 < S0 → post-fix S3 < S1 < S2 < S0

## 2. Paired S3 − baseline (mean ± 95% CI, % change, significant?, S3 better/worse/tied)

### triangle_area

| vs | pre-fix | post-fix | significance change |
|---|---|---|---|
| S0 | -764.7 ± 412.9 (-37.6%, sig, 89/9/2) | -837.2 ± 457.2 (-35.2%, sig, 88/10/2) |  |
| S1 | -255.8 ± 332.7 (-16.8%, n.s., 89/11/0) | -164.3 ± 87.6 (-9.6%, sig, 73/27/0) | **n.s.→sig** |
| S2 | -289.6 ± 148.1 (-18.6%, sig, 70/30/0) | -240.6 ± 140.6 (-13.5%, sig, 69/31/0) |  |

### weighted_triangle_area

| vs | pre-fix | post-fix | significance change |
|---|---|---|---|
| S0 | -1314.3 ± 699.2 (-40.9%, sig, 88/10/2) | -1392.4 ± 748.5 (-38.9%, sig, 90/8/2) |  |
| S1 | -528.5 ± 519.5 (-21.7%, sig, 90/10/0) | -292.5 ± 136.6 (-11.8%, sig, 73/27/0) |  |
| S2 | -510.8 ± 248.4 (-21.2%, sig, 75/25/0) | -412.3 ± 226.1 (-15.9%, sig, 74/26/0) |  |

### gini_restore

| vs | pre-fix | post-fix | significance change |
|---|---|---|---|
| S0 | 0.035 ± 0.008 (+18.8%, sig, 14/80/6) | 0.017 ± 0.005 (+10.5%, sig, 19/71/10) |  |
| S1 | 0.053 ± 0.009 (+31.4%, sig, 10/90/0) | 0.016 ± 0.005 (+10.0%, sig, 27/73/0) |  |
| S2 | 0.033 ± 0.008 (+17.8%, sig, 22/78/0) | 0.017 ± 0.006 (+10.9%, sig, 34/66/0) |  |

### min_time_avg_cri

| vs | pre-fix | post-fix | significance change |
|---|---|---|---|
| S0 | 0.064 ± 0.017 (+22.0%, sig, 78/20/2) | 0.069 ± 0.017 (+23.8%, sig, 81/17/2) |  |
| S1 | 0.094 ± 0.017 (+36.2%, sig, 85/15/0) | 0.022 ± 0.015 (+6.6%, sig, 61/39/0) |  |
| S2 | 0.025 ± 0.015 (+7.7%, sig, 63/37/0) | 0.021 ± 0.015 (+6.2%, sig, 62/38/0) |  |

### p90_access_restore

| vs | pre-fix | post-fix | significance change |
|---|---|---|---|
| S0 | -498.2 ± 275.3 (-29.6%, sig, 80/10/10) | -498.2 ± 275.3 (-29.6%, sig, 80/10/10) |  |
| S1 | -195.4 ± 270.3 (-14.1%, n.s., 88/12/0) | -106.0 ± 54.1 (-8.2%, sig, 76/23/1) | **n.s.→sig** |
| S2 | -190.8 ± 95.5 (-13.9%, sig, 79/21/0) | -134.0 ± 79.2 (-10.1%, sig, 72/28/0) |  |

## 3. Static-score coverage and degeneracy (should be unchanged)

- damaged road links with nonzero S3 score: pre 290/793 (36.6%) → post 290/793 (36.6%)
- S3 road order identical to S0: pre 2/100 → post 2/100

## 4. k-sensitivity (15-scenario batch, S3 only)

| k | pre-fix triangle | post-fix triangle | pre-fix weighted | post-fix weighted |
|---|---|---|---|---|
| 1 | 1695.7 | 2187.7 | 2559.1 | 3126.1 |
| 2 | 1712.4 | 2206.6 | 2608.4 | 3178.0 |
| 3 | 1711.3 | 2205.1 | 2607.0 | 3176.2 |
| 5 | 1676.8 | 2167.3 | 2542.2 | 3107.8 |

post-fix spread across k: 1.8% of mean
