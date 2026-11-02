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

## Stage 2 — HMM bracket
| MATLAB | sqw_ref | notes |
|---|---|---|
| `hmmFitPhased` | `em_fit` | batched Baum-Welch, Dirichlet(1) starts; `Start=` for one fixed start |
| `labelStates`, `isDegenerateH1`, `hmmBic` | `label_states`, `is_degenerate`, `bic` | |
| `hmmViterbiPhased`, `hmmFFBS`, `pairCounts`, `posteriorSummaries` | `viterbi`, `ffbs`, `pair_counts`, `posterior_summaries` | |
| `hmmSimulatePhased` | `simulate` | refuses rows that do not sum to 1 (`fillEmptyRows` first) |
| `roundTripsHidden`, `roundTripStats` | `round_trips_hidden`, `round_trip_stats` | |
| `readStored20bFits`, `check20bExact`, `checkSubset20b` | — | cClaude's stored 20b fits; exact and statistical checks |

Tests: `test_sqw_stage2`.

## Stage 3 — model-free probes (20p, 20q)
| MATLAB | sqw_ref | notes |
|---|---|---|
| `chain2FitPhased`, `pairChainStructure` | `mc2_fit`, `B_PAIR`, `ALLOW_PAIR` | second-order chain as 9 pair states, k = 3(y_{t-1}-1)+y_t |
| `fillEmptyRows` | `sim_rows` | simulation only; also applied to HMM rows never left (sqw_ref jumps to the last state instead) |
| `blockKL`, `blockKLStats` | `block_kl`, `block_kl_stats` | |
| `cycleIncrements`, `turEstimate`, `turExact`, `turExactCT`, `hmmObservedChain` | `D_CYCLE`, `tur_estimate`, `tur_exact(_ct)`, `hmm_observed_chain` | turExact needs stochastic rows; any K in hmmObservedChain |
| `sqwMemoryFly`, `sqwMemoryGroup` | `sqw15.py`, `sqw15b.py`, `sqw15_report.py` | 20p per fly / group: nulls A, B (A2 halves for DD), HMM and MC2 predictions; CSV, figure |
| `sqwTURFly`, `sqwTURGroup` | `sqw16.py`, `sqw16b.py`, `sqw16_report.py` | 20q per fly / group, rule of three (amendment E); CSV, figure |
| `sqwListFiles`, `sqwFlyInput`, `wilcoxonP`, `binomGreater`, `holmAdjust` | — | helpers |

Tests: `test_sqw_stage3`. Raw-data check: `validate_20pq(dataRoot, jsonlFile, outDir)`.

## Figure exports — Plan A Figs A1–A4 (FIGURE_RECIPES section K)
Each export writes the plotted per-fly points and the per-block medians / CIs as CSV (for Origin) and .mat, and
draws the figure (theme pinned light, genotype order CS → CS DD → per0 LD → per0 DD → fumin → sss, light before dark)
as PNG and vector PDF.

| MATLAB | Python | notes |
|---|---|---|
| `planAData` | `sqw17.py` | six `sqwObservableGroup` outputs in paper order (sss = the 54-file 20s-A set) plus the HMM triangles; stored 20b fits matched to files BY NAME |
| `planAFlyHMM` | `sqw17.py` "hid", `sqw2.decode` | HMM count J and EPR, mean over 200 FFBS paths, per fly x phase |
| `planALikelihoods` | `sqw13.py`, `sqw14.py` | per-fly 20n / 20o log-likelihoods (MC1, HMM, MC2 at 12 h; MC1, HMM at 6 h), dBIC, R_mem, R_drift |
| `exportFigA1` .. `exportFigA3` | `fig_planA.py` | direction and flux (main text observed-level only; `HmmMarkers=true` for the SI version); EPR vs J and Onsager L; energy landscape |
| `routeDurationObserved`, `routeDurationD` | `route_duration_observed`, `route_duration_D` | 20k route-vs-duration records and D (golden G11) |
| `exportFigC7` | `fig_R1_obs.py`, `sqw10_report.py`, `route_bins.py` | Fig C7, sleep-length memory, observed only: recomputes 20k observed R1 from raw data (stars from the ledger by default) |
| `exportFigA4` | `fig_beneath.py` | reads the per-fly CSVs `sqwMemoryGroup` / `sqwTURGroup` write, plus `planALikelihoods` |
| `exportFigC10` | `fig_clockface.py` | Fig C10 clock face, from `planB_20wB.json` (20w-B); the cosine is refitted, its amplitude checked against the stored AJ |
| `exportFigC12` | `fig_kir.py` | Fig C12 Kir2.1 silencing, from `y_kir.json` (20y); ledger Holm stars copied, not recomputed |
| `planAStyle`, `planAFigure`, `planAStrip` | — | colours, fonts, blocks; figure with theme pinned; strip panel |

Tests: `test_sqw_figA`.

Fig C8 (per0 DD and CS DD wake hazards at 12/16/20/24 bins) is a bout-suite figure, not an sqw one:
`analysis_20r/figC8_wake_hazard.m`, run by name from the Icode folder; it recomputes the ledgered rise tests from the
saved 20r bouts and checks every RiseNullP against the ledger. Outputs in `analysis_20r/figC8`.

Added 2026-11-02 (FIGURE_RECIPES notes dated 2026-11-02):

| MATLAB | Python | notes |
|---|---|---|
| `exportFigA1` / `exportFigA2` | `fig_planA.py` | sss mean marker in A1 (a, b); A2 sss title note and panel (c) A_eff = EPR/\|J\| |
| `exportFigSssReversal` | `fig_sss_reversal.py` | per-fly chain J beside HMM J, CS and sss; no-cycle flies (|J| < 1e-12) plotted at 0, left out of the % and the test |
| `exportFig1Emissions` | `fig_emission.py` | median labelled B3 per genotype |
| `exportFigDwellDirection` | `fig_dwell_direction.py` | CS night, direction-dependent Q and W dwells |
| `exportFigC13Profiles` | `fig_planB_21.py` | C13 / PRESPEC 21 profiles; checks itself against fig_planB21_profiles.json |
| `routeDurationHidden`, `planARouteHidden`, `exportFigC7SI` | `route_duration_hidden`, `sqw10.py`, `fig_R1_main.py` | SI memory figure with the hidden panel; the hidden 20k null takes ~1 h on 2 workers (cache it) |
