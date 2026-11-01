# sqw — MATLAB port of the S/Q/W analysis (Plan A)

MATLAB (R2025a+, tested in R2026a) port of `hmm_audit/sqw_ref.py`, the only code to port
(`hmm_audit/START_HERE.md`, `PORTING_GUIDE.md`). Branch `matlab-port`. Each function names the
`sqw_ref` function it ports in its first help line.

## Conventions (differ from the Python and the golden files)
- **1-based everything.** Symbols `1 = Q`, `2 = W`, `3 = S`, `0 = missing`, `4 = wake of unknown type`
  (Python 0, 1, 2, -1, -2). Phases `1..nP` (light 1, dark 2; 6-h blocks 1..4). Bout starts and minutes 1-based.
- **Transition arrays are `K x K x nP`** (`P(i,j,p)`); `jsondecode` of a golden file gives `nP x K x K` —
  convert with `permute(P, [2 3 1])`.
- **Phase of arrival.** `ph(t)` is the phase of the step INTO minute `t`; every count goes to its arrival minute.
- **Sign.** Forward cycle W -> Q -> S -> W positive (count J, chain J, A, c_obs, c_pred).
- **No toolbox in the core.** t quantiles via `betaincinv` (`tQuantile`), percentiles by numpy's rule
  (`quantile7`; MATLAB's `prctile` uses a different rule), Spearman with average ranks (`spearmanRho`).
- **Random numbers** come from a private `RandStream` (`RandomSeed` options), never the global `rng`.
  numpy and MATLAB streams differ, so bootstrap intervals agree with the Python only within Monte Carlo error.

## Reused, not rewritten (byte-identical copies from `Box\MATLAB scripts\Rockefeller Matlab scripts`)
| file | SHA-256 | use |
|---|---|---|
| `epr_count_based.m` | `5836d22d…52e90c9` | count EPR, rule-of-three bound (`'bound'`); takes `nSteps` in SECONDS: pass `hours*3600` |
| `stationary_dist.m` | `2222e10a…ee8b3` | stationary distribution by the QR solve PORTING_GUIDE prefers |
If either is changed in the Rockefeller folder, re-copy it here and re-run the tests.

## Stage 1 — observable core
| MATLAB | sqw_ref | notes |
|---|---|---|
| `readCounts` | `read_counts` | refuses headers, negatives, fractions |
| `sleepBouts` | `sleep_bouts` | `[start len]`, 1-based; lab bridging (b) |
| `wakeBoutEvents`, `wakeReference` | `wake_bout_events`, `wake_reference` | struct array of events; r_f per phase |
| `sqwSymbols` | `sqw_symbols` | `UnknownWake="missing"` (ledger) or `"unknown"` (recommended, code 4) |
| `minutePhases` | `phases` | `"12h"`, `"6h"`, `"none"` |
| `hmmEmission`, `hmmForward`, `hmmEmStep` | `emission_lik`, `forward`, `em_step` | needed by the chain fit; G3 / G9 exact |
| `chainFitPhased` | `mc_fit` | observed first-order chain, missing minutes marginalised |
| `observablePairCounts`, `cObsMinutes`, `jCycleCount` | same names | |
| `modelEnergetics`, `routingPredicted`, `energyLandscape` | `model_energetics`, `xn/yn/c_pred`, `landscape` | Inf kept for one-way edges |
| `routeEvents`, `routeCounts`, `routeCountsObserved`, `xyCounts` | `route_events`, `route_counts`, `routes_observed`, `xy_counts` | |
| `observableQuantities` | `observable_quantities` | per phase; adds hours, N, routes |
| `onsagerFit`, `onsagerBootstrap`, `bootstrapDifferenceP`, `eprJModels` | same | |
| `sqwObservableFly`, `sqwObservableGroup`, `sqwOnsagerGroup` | `sqw17.py`, `sqw17_report.py` | drivers: structs, CSV for Origin, figures |
| `plotSqwFly`, `plotSqwGroup`, `plotSqwOnsager` | — | figures (theme pinned light) |
| `medianCI`, `quantile7`, `tQuantile`, `spearmanRho` | — | helpers |

Tests: `test_sqw_stage1` (run by name from this folder). Raw-data check: `validate_20sA(dataRoot, outDir)`.
Validation record: `hmm_audit/PORT_VALIDATION.md`.
