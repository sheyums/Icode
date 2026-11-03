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
| `readStoredHmmFits`, `checkHmmLikelihoodExact`, `checkHmmRefitSubset` | — | cClaude's stored 20b fits; exact and statistical checks |

Tests: `test_sqw_hmm`.

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

Tests: `test_sqw_modelFree`. Raw-data check: `validateMemoryTURLedger(dataRoot, jsonlFile, outDir)`.

## Figure exports — Plan A Figs A1–A4 (FIGURE_RECIPES section K)
Each export writes the plotted per-fly points and the per-block medians / CIs as CSV (for Origin) and .mat, and
draws the figure (theme pinned light, genotype order CS → CS DD → per0 LD → per0 DD → fumin → sss, light before dark)
as PNG and vector PDF.

| MATLAB | Python | notes |
|---|---|---|
| `buildPerFlyData` | `sqw17.py` | six `sqwObservableGroup` outputs in paper order (sss = the 54-file 20s-A set) plus the HMM triangles; stored 20b fits matched to files BY NAME |
| `hmmPosteriorFluxEPR` | `sqw17.py` "hid", `sqw2.decode` | HMM count J and EPR, mean over 200 FFBS paths, per fly x phase |
| `chainVsHmmLikelihoods` | `sqw13.py`, `sqw14.py` | per-fly 20n / 20o log-likelihoods (MC1, HMM, MC2 at 12 h; MC1, HMM at 6 h), dBIC, R_mem, R_drift |
| `exportDirectionFluxFigure` .. `exportLandscapeFigure` | `fig_planA.py` | direction and flux (main text observed-level only; `HmmMarkers=true` for the SI version); EPR vs J and Onsager L; energy landscape |
| `routeDurationObserved`, `routeDurationD` | `route_duration_observed`, `route_duration_D` | 20k route-vs-duration records and D (golden G11) |
| `exportSleepMemoryFigure` | `fig_R1_obs.py`, `sqw10_report.py`, `route_bins.py` | Fig C7, sleep-length memory, observed only: recomputes 20k observed R1 from raw data (stars from the ledger by default) |
| `exportHiddenVsObservedFigure` | `fig_beneath.py` | reads the per-fly CSVs `sqwMemoryGroup` / `sqwTURGroup` write, plus `chainVsHmmLikelihoods` |
| `exportClockFaceFigure` | `fig_clockface.py` | Fig C10 clock face, from `planB_20wB.json` (20w-B); the cosine is refitted, its amplitude checked against the stored AJ |
| `exportKirSilencingFigure` | `fig_kir.py` | Fig C12 Kir2.1 silencing, from `y_kir.json` (20y); ledger Holm stars copied, not recomputed |
| `paperFigureStyle`, `newPaperFigure`, `plotStripPanel` | — | colours, fonts, blocks; figure with theme pinned; strip panel |

Tests: `test_sqw_figures`.

Fig C8 (per0 DD and CS DD wake hazards at 12/16/20/24 bins) is a bout-suite figure, not an sqw one:
`analysis_20r/wakeHazardFigure.m`, run by name from the Icode folder; it recomputes the ledgered rise tests from the
saved 20r bouts and checks every RiseNullP against the ledger. Outputs in `analysis_20r/figC8`.

Added 2026-11-02 (FIGURE_RECIPES notes dated 2026-11-02):

| MATLAB | Python | notes |
|---|---|---|
| `exportDirectionFluxFigure` / `exportEprOnsagerFigure` | `fig_planA.py` | sss mean marker in A1 (a, b); A2 sss title note and panel (c) A_eff = EPR/\|J\| |
| `exportSssReversalFigure` | `fig_sss_reversal.py` | per-fly chain J beside HMM J, CS and sss; no-cycle flies (\|J\| < 1e-12) plotted at 0, left out of the % and the test |
| `exportHmmEmissionsFigure` | `fig_emission.py` | median labelled B3 per genotype |
| `exportDwellDirectionFigure` | `fig_dwell_direction.py` | CS night, direction-dependent Q and W dwells |
| `exportWakeEndProfilesFigure` | `fig_planB_21.py` | C13 / PRESPEC 21 profiles; checks itself against fig_planB21_profiles.json |
| `routeDurationHidden`, `routeDurationHiddenNull`, `exportSleepMemoryHiddenFigure` | `route_duration_hidden`, `sqw10.py`, `fig_R1_main.py` | SI memory figure with the hidden panel; the hidden 20k null takes ~1 h on 2 workers (cache it) |

## Batch tools for new data — lab coding (added 2026-11-02)
Every file these write uses the lab's coding **S = 0, Q = 1, W = 2, missing = NaN** (sleep stays 0, as in
`sleeptrace`). Internally the port keeps 1 = Q, 2 = W, 3 = S, 0 = missing, which MATLAB needs to index by symbol;
`sqwToLabCode` / `labCodeToSqw` convert, and `sqwObservableFly` returns both (`F.symbols`, `F.y`).

| MATLAB | does | notes |
|---|---|---|
| `sqwBatchSymbols(folder)` | count files -> one `<fly>_symbols.csv` per fly (minute, zt_h, phase, symbol, active_frac, sleep_s) + `symbols_summary.csv` | `LightDark=false` for DD; theta 0.5; output folder `<folder>_symbols` BESIDE the input. Scoring = `sqwSymbols`, NOT `sleeptrace` + `ironout` |
| `sqwBatchHMM(folder)` | K = 2, 3, 4 per fly on symbol files (or count files): ll, BIC (n = minutes, n = changes), selected K, PRESPEC 20 flag, labelled K = 3 B and P, dominant symbol per state for K = 2 / 4 | `<Name>_perfly.csv`, `<Name>.mat`, `<Name>_order.png`; per-fly fits cached in `<Name>_fits/` so an interrupted batch resumes; ~3 min per fly per worker (20 starts) |
| `plotSqwHMMOrder(R)` / `({R1,...}, names)` | K = 3 emissions, BIC selection (both n), K = 4 "states emitting mostly S" | one row for one batch, grid for several |

Tests: `test_sqw_batch` (12). Validation against the stored 20b fits: `hmm_audit/PORT_VALIDATION.md`.

Typical use:
```matlab
X = sqwBatchSymbols("D:\data\CS");                 % writes D:\data\CS_symbols\...
R = sqwBatchHMM(X.outDir);                          % writes hmm_perfly.csv, hmm.mat, hmm_order.png there
Xdd = sqwBatchSymbols("D:\data\per0DD", LightDark=false); Rdd = sqwBatchHMM(Xdd.outDir);
plotSqwHMMOrder({R, Rdd}, ["CS" "\itper^0\rm DD"]);
```

## Renamed 2026-11-02 (meaningful names; old names appear in PORT_VALIDATION.md and FIGURE_RECIPES.md)
Output file names (figA1_*.png, figC7.mat, ...) are unchanged.

| old | new | what it is |
|---|---|---|
| `exportFigC13Profiles` | `exportWakeEndProfilesFigure` | C13 / PRESPEC 21 profiles |
| `exportFigC7SI` | `exportSleepMemoryHiddenFigure` | Fig C7 SI (hidden level) |
| `exportFigC7` | `exportSleepMemoryFigure` | Fig C7 |
| `exportFigC10` | `exportClockFaceFigure` | Fig C10 |
| `exportFigC12` | `exportKirSilencingFigure` | Fig C12 |
| `exportFig1Emissions` | `exportHmmEmissionsFigure` | Fig 1 HMM emission panel |
| `exportFigDwellDirection` | `exportDwellDirectionFigure` | dwell-direction figure |
| `exportFigSssReversal` | `exportSssReversalFigure` | sss reversal figure |
| `exportFigA1` | `exportDirectionFluxFigure` | Fig A1 |
| `exportFigA2` | `exportEprOnsagerFigure` | Fig A2 |
| `exportFigA3` | `exportLandscapeFigure` | Fig A3 |
| `exportFigA4` | `exportHiddenVsObservedFigure` | Fig A4 |
| `planAData` | `buildPerFlyData` | per fly x phase data behind Figs A1-A3 (was Plan A data) |
| `planAFigure` | `newPaperFigure` | figure with the paper style, theme pinned light |
| `planAStyle` | `paperFigureStyle` | colours, order, font |
| `planAStrip` | `plotStripPanel` | one strip panel |
| `planAFlyHMM` | `hmmPosteriorFluxEPR` | HMM flux and EPR over posterior paths |
| `planALikelihoods` | `chainVsHmmLikelihoods` | 20n/20o per-fly likelihoods |
| `planARouteHidden` | `routeDurationHiddenNull` | hidden-level 20k null |
| `check20bExact` | `checkHmmLikelihoodExact` | 20b exact likelihood check |
| `checkSubset20b` | `checkHmmRefitSubset` | 20b refit of a subset |
| `readStored20bFits` | `readStoredHmmFits` | cClaude stored 20b fits |
| `validate_20pq` | `validateMemoryTURLedger` | 20p/20q ledger check |
| `test_sqw_figA` | `test_sqw_figures` | figure-export tests |
| `test_sqw_stage2` | `test_sqw_hmm` | HMM tests (was stage 2) |
| `test_sqw_stage3` | `test_sqw_modelFree` | model-free probe tests (was stage 3) |
| `figC8_wake_hazard` | `wakeHazardFigure` | Fig C8 (analysis_20r) |
