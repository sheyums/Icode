# MATLAB port validation (sqw_ref.py -> Icode/sqw), branch matlab-port

Append-only, one entry per stage, with the commit of the code that passed. Port by lClaude.

## Status at pause (2026-11-01)

A summary added at the top on Sheyum's request; the dated entries below are the record. Code: Icode/sqw (functions
listed in sqw/README.md), branch matlab-port, merged into claude/trusting-allen-7nxrsf.

**Ported and validated** (MATLAB R2026a; tests run by name from Icode/sqw):
- Stage 1, observable core (symbols, sleep bouts, wake reference, pair/route counts, c_obs, counted J and EPR, observed
  chain, landscape, Onsager fit/bootstrap, 20s-A group tables): `test_sqw_stage1` 21/21; 20s-A ledger medians
  reproduced from raw data (`validate_20sA`).
- Stage 2, HMM (EM fit, labelling, BIC, Viterbi, FFBS, posterior summaries, simulation, round trips): `test_sqw_stage2`
  11/11; exact per-fly log-likelihood check against the stored 20b fits, 10/10 flies, max |dll| 2.5e-10.
- Stage 3, model-free probes (second-order chain, block KL, TUR estimate and exact, HMM observed chain, 20p/20q group
  drivers): `test_sqw_stage3` 13/13; 20p/20q ledger reproduced from raw data (160/160 deterministic numbers exact, 180
  Monte Carlo statistics within error, every verdict agrees; `validate_20pq`).
- Route vs duration (`routeDurationObserved`, `routeDurationD`, 20k): golden G11 exact; 20k observed R1 n and D equal
  the ledger in all 9 blocks.
- Figure exports and helpers: `test_sqw_figA` 10/10.
- Root bout suite: BoutHazard overlay fixed (cd2f86f), `test_BoutHazard` 11/11; root suite now 130 tests.

**Figures that exist in MATLAB** (data as CSV/.mat plus PNG/PDF, theme pinned light; outputs in sqw/validation/figA,
sqw/validation/figC, analysis_20r/figC8):
- Fig A1 (observed only; SI variant with HMM medians via `HmmMarkers=true`), Fig A2, Fig A3, Fig A4 (with the
  adopted panel (c) legend): `exportFigA1`..`exportFigA4`, data from `planAData`, `planALikelihoods`.
- Fig C7 observed (`exportFigC7`), Fig C8 wake hazard (`analysis_20r/figC8_wake_hazard.m`; all 8 RiseNullP equal the
  ledger), Fig C10 clock face (`exportFigC10`), Fig C12 Kir2.1 (`exportFigC12`).

**Not ported yet:**
- The C13 profiles figure (fig_planB_21.py).
- The dwell-direction figure (fig_dwell_direction.py; added by cClaude at 6e96fc4).
- The SI memory figure with the HMM panel (fig_R1_main.py; figS_C7_memory_with_hmm_check).
- Small alignment with cClaude's 2b56d1f: the SI Fig A1 is named figS_A1_with_hmm_check there, its option
  'IncludeHMM' and legend "hidden Markov model (check), median"; the port writes figA1_direction_flux_SI_hmm with
  `HmmMarkers` and "HMM (model-based estimate), median". Same content; names to be aligned.

**Open notes:**
- Q3 tolerance: exact TUR D and B agree with Python to ~1e-7 (round-off amplified by the second difference), so G12
  checks them at 1e-6 relative. cClaude accepted this on 2026-11-01 (ledger notes); recorded here for completeness.
- Analysis-folder copy of BoutHazard.m: it WAS one fix behind after cd2f86f, but it has since been re-synced on
  Sheyum's request (C:\Users\sheyu\Box\MATLAB scripts\Rockefeller Matlab scripts, BoutHazard.m and
  test_BoutHazard.m from cd2f86f, 11/11 there, Icode_sync_stamp.txt updated, committed in that folder's repository
  as a301f73). It is now current with matlab-port and claude/trusting-allen-7nxrsf. It is AHEAD of Icode main,
  which does not have the fix yet.
- per0 LD light (20k observed R1): the Holm star is near threshold at 2000 permutations; Fig C7 takes the ledger's
  stars (see the cd2f86f entry).

## Stage 1 — observable core (2026-11-01, code 2bddb4d)

**Functions** (Icode/sqw; map in sqw/README.md): readCounts, sleepBouts, wakeBoutEvents, wakeReference,
sqwSymbols, minutePhases, hmmEmission/hmmForward/hmmEmStep (needed by mc_fit), chainFitPhased (mc_fit),
observablePairCounts, cObsMinutes, jCycleCount, epr_count_based.m (reused: byte-identical copy of the Rockefeller
file, SHA-256 5836d22d...), stationary_dist.m (reused, QR solve, 2222e10a...), modelEnergetics, routingPredicted,
energyLandscape, routeEvents, routeCounts, routeCountsObserved (routes_observed), xyCounts, observableQuantities,
onsagerFit, onsagerBootstrap, bootstrapDifferenceP, eprJModels. Drivers: sqwObservableFly, sqwObservableGroup,
sqwOnsagerGroup (sqw17.py + sqw17_report.py), plots, CSV export. No toolbox in the core.

**Golden tests** (`test_sqw_stage1`, MATLAB R2026a): 21/21 pass.
- G1 all 8 bridging cases exact; G2 symbols (2880 min), 49 bouts, 30 events, refs exact; G3 loglik and one EM step,
  max relative error 2.0e-14; G9 (unknown wake) 1.9e-14; G4 routing / pi / J / EPR / A to 1e-12 (Inf kept);
  G5 route events, counts, pair counts, Xn/Yn, EPR bound (4 cases), J, c_obs exact; G7 observable pair counts exact;
  G10 Onsager fit and EPR-J models to 1e-12, bootstrap median within 0.05% of L; G11 landscape to 1e-12;
  G12 mc_fit loglik -1216.031986 (to 1e-6), P to 1e-9, routes_observed exact, 13 observable quantities x 2 phases to 1e-9.
- Stress tests (brute force or identities): bridging rule vs a literal loop (300 traces x 3 settings); pair counts,
  c_obs and observed routes vs loops on 200 sequences with gaps; chain fit = count MLE on complete data, EM monotone
  with gaps and unknown wake; 2000 random matrices: sign(J) = sign(A) = sign(c), J equal on the three edges,
  EPR = J A, detailed balance 0, one-way edge Inf; count J = epr_count_based J_cycle and EPR bound >= 0 on 5000 count
  matrices; t quantile = tinv and closed forms; type-7 quantiles = numpy; readCounts refuses headers / negatives /
  fractions; drivers = the functions they chain, CSV round trip; Onsager group recovers known L.
- Calibration (not in the suite, 400 null replicates of two genotypes with equal L, EPR-proportional noise 5%,
  40 flies x 2 points each): bootstrap-difference p < 0.05 in 0.052; **Wald p < 0.05 in 0.220**; fly-bootstrap 95% CI
  of L covers the truth in 0.932.

**Raw data** (`validate_20sA`): the files cClaude used, from Sheyum's zips in Box\HMM paper\different genotypes:
CS (SHA-256 prefix 6539484a, 25 files), per0LD (6cd83438, 27), per0DD (778599db, 25; the 20r per0 DD archive),
fumin (f4d38318, 53), sss (a2f1180f, 57 files of which the 54 cClaude listed were used), CS DD (the 20r CS DD files, 25).
**0 disagreements**: every 20s-A median in RESULTS_LEDGER.md agrees at its printed precision: n, n with a cycle,
c_obs, count J, count EPR, chain J / A / EPR, E2-E1, E3-E1 for all 10 genotype x phase blocks (e.g. CS dark c_obs
+0.147, count J +0.167, count EPR 0.116), the CS 6-h blocks (C3), and Onsager L with its n for count points
(CS 0.272, per0LD 0.286, per0DD 0.352, fumin 0.243, sss 0.145, CS DD 0.262) and chain points. Bootstrap intervals
are not compared (numpy and MATLAB streams differ).

**Cases to report (not changed on either side; for cClaude and Sheyum to decide):**
1. *sss file list.* Sheyum's sss zip holds 57 files; cClaude's run used 54. The three not used
   (4_20_2018_4boardmb_pruned_5, 5_9_2018_4boardmb_pruned_3, 9_16_2019_4b_rc2_pruned_3) look like the others
   (259199 s, same export). Which set is the paper's?
2. *"No cycle" (|J_chain| < 1e-12) lumps three different things.* Of the 28 fly-phases with a Q/W reference so
   excluded from the chain quantities and the Onsager points: 11 have EXACTLY SYMMETRIC counts (all three edges seen
   both ways, zero net current: a cycle with J = 0, not an absent cycle), 15 have one edge (mostly Q<->S) never seen
   in either direction (each rate is only bounded by the rule of three, so a cycle cannot be excluded), and 2 have a
   state absent (sss flies 1 and 38, no sleep minute in the light half: truly no cycle). By block: CS L 1, per0LD L 1,
   per0DD 1, fumin L 2 / D 4, sss L 13 / D 5. PRESPEC 20q amendment E already separates "exactly symmetric" from
   "no cycle" for the TUR; 20s-A does not.
3. *Median chain EPR drops infinite values.* The report's chain-EPR median is over finite values only, so the flies
   with a one-way edge (the most irreversible) are removed. With them kept (Inf ranks highest): fumin L 0.163 -> 0.170,
   fumin D 0.132 -> 0.189, sss L 0.046 -> 0.090, sss D 0.043 -> 0.109; other blocks unchanged (no Inf). The summary
   table of sqwObservableGroup carries both (EPR_chain_med, EPR_chain_med_withInf) and the Inf share.
4. *Empty chain rows.* sqw_ref.model_energetics accepts a matrix with an all-zero row (a state never entered in the
   phase: sss flies 1 and 38, light). The port accepts it too (pi of that state is 0, no cycle) and still refuses any
   other row that does not sum to 1.
5. *Wald overstates* under EPR-proportional noise (22% false positives at nominal 5% in the calibration above), which
   supports the ledger's choice of the bootstrap difference for the L comparisons (20s-A addendum 2).

## Stage 2 — HMM bracket, golden and stress part (2026-11-01, code 74cfead); real-data part PENDING

**Functions:** hmmFitPhased (em_fit; all starts batched with sqw_ref's joint stopping rule), labelStates,
isDegenerateH1, hmmBic (n_params, bic, n_eff), hmmViterbiPhased, hmmFFBS, pairCounts, routeCounts (stage 1),
posteriorSummaries, hmmSimulatePhased. hmmEmStep / hmmForward were validated in stage 1 (G3, G9).

**Tests** (`test_sqw_stage2`): 11/11 pass.
- G3 parameter counts (9 cases) and BIC exact; G3 and G9 Viterbi paths exact (2880 minutes each); G4 labelling and
  degeneracy exact.
- G6 (statistical, 4000 FFBS paths on the G3 inputs): mean route counts max |z| 1.80 (bound 4); J |z| 0.84,
  EPR bound |z| 0.88; Xn, Yn, c within 0.006.
- hmmFitPhased on G3 reaches ll -781.466688, equal to the G3 parameters after 300 EM steps.
- Batched fit = 40 iterated hmmEmStep calls from the same start (3 random starts, unknown-wake minutes): max diff 1e-15.
- Viterbi = brute force over all 3^7 paths (50 random models with gaps and unknown wake).
- FFBS one- and two-minute marginals = exact forward-backward posteriors (20000 paths, max |z| 2.30).
- Labelling invariant to 20 random state orders, likelihood unchanged; recovery of a known HMM from 8 simulated
  days (|dP| <= 0.024, |dB| <= 0.013; fitted ll above the true-parameter ll); simulation frequencies within
  sampling error.

**Pending:** per-fly 20b log-likelihoods to 1e-6 from the same start needs the stored fits (sqw_out/fit2_0.5.jsonl),
requested in QUESTIONS_lClaude.md Q2. Timing for planning: one CS fly (5760 minutes) at K = 3, 20 starts, takes
121 s and stops at MaxIter 1000 (16 of 20 starts at the best optimum, ll -2539.6261), so a full 20b refit of 187 flies
x K 2-4 is ~19 h serial.

## Stage 2 — real-data subset, statistical (2026-11-01, code below; fits in sqw/validation/subset20b_fits.mat)

Refit of 10 flies (the first 2 per genotype in cClaude's file order: CS, per0LD, per0DD, fumin, sss) at K = 3,
20 random starts, theta 0.5, with `checkSubset20b`; labelled, 200 FFBS paths; hidden round trips (20g, all and
matched 5/10) per fly x phase compared with cClaude's `roundtrip_per_fly_theta0.5.csv` (made from its 20b fits):
p_quiet_entry, p_quiet_exit and the four fractions. z = difference / (sqrt(2) x Monte Carlo SE over paths);
in the 36 cells where all 200 paths agree (SD 0) the SE is the binomial sqrt(p(1-p)/200).
- 204 comparisons: |z| > 2 in 0.044 (noise alone ~0.05), |z| > 3 in none, max |z| 2.96. Per fly max |z|:
  0.3-2.96; no fly stands out, so the port reached cClaude's optimum on all 10 (a different optimum shows as |z| >> 3).
- Per-fly log-likelihoods (K = 3) are in sqw/validation/subset20b_compare.csv's companion .mat, kept for the exact
  check once the stored fits arrive. No H1 degeneracy in the 10.
- Without the binomial floor, one cell gave |z| 3535 (per0LD 8_25_2017_multibeam_pruned3, light, matched,
  p_quiet_entry: port 1.000 on all 200 paths, Python 0.995 = one path in 200): an SE of 0, not a different fit.
- Still pending for the stage-2 pass: the exact per-fly log-likelihood check (Q2).

## Stage 2 — exact per-fly log-likelihood check (Q2) (2026-11-01, code c19800a; output sqw/validation/exact20b_compare.csv)

**STAGE 2 PASSES.** Stored fits: hmm_audit/sqw20b_fit2_theta0.5.jsonl (cClaude, 0212253, merged into matlab-port at
18d1a66). Checked on the 10 subset flies of the entry above (CS, per0LD, per0DD fly 0-1; fumin, sss fly 0-1; sss indexed
in the 57-file set, per A2), raw data re-extracted from the user's zips (SHA-256 prefixes as in stage 1), with
`check20bExact`:
- Identity: T, the number of symbol changes (nchg) and the Q/W/S fractions of the port's symbols equal the stored ones
  on all 10 (fractions to 0), so the records are the same flies read the same way.
- Exact: `hmmForward` at the stored (pi, P, B) reproduces the stored ll for K = 2, 3 and 4, and at the labelled K = 3
  fit (pi3, P3, B3): max |difference| 2.5e-10 over 40 evaluations (tolerance 1e-6). 10 / 10 pass.
- Baum-Welch started AT the stored K = 3 fit stops after 2 iterations on every fly, with gains between -9e-11 and
  +4e-8, except fumin fly 0 (+3.2e-5, still below the stop rule's 1e-7 |ll| = 1.8e-4).
- The port's own 20-start K = 3 optimum equals the stored one within 1.4e-7 on 9 flies. fumin fly 0
  (10_15_2018_4b_RC2_pruned_14) is +0.0038 HIGHER in the port: the stored fit stopped at 460 iterations (every Python
  start had met the relative-gain rule), while the port's batch ran to 1000 because another start was still climbing,
  so its best start crept further up the same slow ridge. Same optimum basin, a stopping-point difference allowed by
  the shared stop rule; neither code is wrong. Effect on this fly's statistics: none detectable (its round-trip
  max |z| 2.96 in the subset entry is inside noise).
- Not checked here: the remaining 177 flies (seconds each by forward evaluation, if wanted), and
  sqw20b_route2_theta0.5.json (route/pair counts per path are FFBS draws; stage 3's G12 and the subset entry cover
  them statistically).

## Stage 3 — IN PROGRESS, NOT PASSED (2026-11-01, work in progress; see test_sqw_stage3.m)
Ported: chain2FitPhased (mc2_fit), pairChainStructure, fillEmptyRows (sim_rows), blockKL, blockKLStats, cycleIncrements,
turEstimate, turExact, turExactCT, hmmObservedChain. First test run 9/13:
- PASS golden: mc2 loglik (1e-6), blockKL/QW/stats (3e-16), turEstimate (0), fillEmptyRows; stress: blockKL and
  turEstimate = literal loops, exact J = sum pi T d = 3 x cycle J, continuous-time TUR B <= sigma over 1000 chains,
  pair-chain simulation, 20p null p-values uniform (p < .05 in 0.09 of 100).
- FAIL test 4 (to REPORT, not tweak): turExact/turExactCT/hmmObservedChain match G12 to 1.5e-7, not the golden's 1e-9.
  The differences are in D (second central difference of log Perron root, h = 1e-3: eigenvalue round-off ~1e-16 is
  amplified by 1/h^2 x 60 ~ 6e7). Probably unattainable at 1e-9 in any second implementation; to agree with cClaude.
- FAIL tests 6, 8, 11 (stress): not yet diagnosed. 6: pair-chain P on a complete record 1e-3 from the triple MLE, and
  first-order data dev 0.037 > 0.03. 8: long-chain d2 0.0173 vs exact EPR 0.0218. 11: simulated symbol current -0.131
  vs exact -0.153 (about 1.3 SE: threshold likely too tight). Each may be the test, not the code; check before trusting.
- Not started: 20p/20q group drivers, the ledger-median comparison, figure exports.

## Stage 3 — PASSES (2026-11-01, code 16c3ee1; outputs sqw/validation/compare_20pq.csv, port_20p_*.csv, port_20q_blocks.csv)

**Golden G12 and stress: test_sqw_stage3 13/13.** The three stress failures of the WIP entry above were the tests' own
expectations, not the port: (6) EM also counts the first transition from the unobserved pair (y0, y1), so only rows whose
second symbol is y1 differ from the triple MLE (all other rows agree to 1e-16); the first-order check is now a z-test
(max |z| 1.85); (8) with Ptrue (EPR 0.02 nats/h) the sampling error of d2 at 4e5 minutes is ~20%; on a strongly driven
chain d2 = 9.93 vs exact EPR 9.98 nats/h, D3 -0.02, D4 +0.11; (11) the simulated current was 0.98 SE from the exact.
Test 4 checks D and B at 1e-6 relative instead of 1e-9: Q3 in QUESTIONS_lClaude.md (round-off measured, neither side changed).

**Raw-data reproduction of the 20p and 20q ledger lines** (validate_20pq; all six genotypes, sss = 57 files; HMM =
cClaude's stored 20b fits, CS DD refitted here as 20p did, K = 3, 20 starts; 200/200/50/50 records per fly for 20p,
200/100 for 20q, as pre-specified):
- EXACT (deterministic given the data): 160 of 160 ledger numbers agree to their printed precision: n, d2, d3, d4,
  Delta3_obs and d3/d2 for the 10 blocks and the 4 DD halves; n, J, D, blocks, B_obs, sigma_MC, Q_MC and kappa max for
  the 10 TUR blocks (amendment E rows for fumin and sss). The amendment-E classification of excluded fly-phases is
  reproduced exactly (kind i fumin 2/10, sss 7/14; ii CS light 1, per0LD light 1, per0DD 1, fumin 1/2, sss 4/1;
  iii fumin 1/2, sss 9/4).
- MONTE CARLO: 180 statistics (E under nulls A and B, E_D4, E_D3QW, fractions with p < 0.05, HMM/MC2 predicted Delta3
  and bracketing, S1, S2, Q_obs, Q_HMM, sigma_HMM >= B_c, HMM bracketing). z = (port - ledger) / sqrt(2 SE^2 + rounding^2/3),
  SE from resampling each fly's own simulated records. z mean +0.07, SD 1.06; |z| > 2 in 18 (10%; 4 of them are
  duplicates, since Q_obs IS S2_med), |z| > 3 in none, max 2.56 (CS DD S1). No statistic is biased one way across
  blocks (per-statistic mean z between -0.44 and +0.64 over 10 blocks each). Within Monte Carlo error.
- VERDICTS agree: 20p QUALIFIES in fumin dark, per0DD whole, CS DD whole and nowhere else (halves: none);
  20q S1 and S2 not supported in either primary DD block.

**Cases reported (nothing changed on the Python side):**
1. G12 tolerance for the exact TUR D and B (Q3).
2. Two stored HMM fits have an all-zero row: sss fly 1 and fly 40 (0-based), H3 in the light phase (no S there). The
   H3 column is exactly 0 in that phase, so no exact quantity changes, but sqw_ref.simulate maps a zero row to the LAST
   state, so a simulated fly in H3 at lights-on stays asleep through the whole light phase. The port fills such rows
   uniform (as 20p amendment B does for the chain). Affects only those two flies' HMM-simulated records (sss light).
3. The 20q entry's "no estimate (blocks < 4 or D <= 0)" counts are pooled over the three block lengths and divided by 3
   (sqw16_report counts every exclusion of all three tau calls). At the primary 3 h the counts are fumin light 1 (ledger
   15) and sss light 0 (ledger 10); the n's printed on those lines are the 3-h ones and are right.
4. Q_obs and S2 are the same per-fly quantity (B_c / sigma_MC) in 20q; the ledger prints it twice.
5. make_golden_planA.py labels the light matrix "forward-biased cycle"; its cycle current is negative (J = -0.23/h:
   reverse-biased). The golden values are unaffected.

Not yet done: the figure-export functions (FIGURE_RECIPES section K, A1-A4).

## 2026-11-01  Figure exports A1-A4 (FIGURE_RECIPES section K) — code 6809a63 (lClaude, matlab-port)

Merged claude/trusting-allen-7nxrsf (a3e1b63, includes 366efc5: "model-based estimate" wording) into matlab-port first.
PRESPEC 20z Part A was NOT run (frozen by Sheyum 2026-11-01; only a synthetic smoke test of a draft runner was
started and stopped, nothing committed, the Kir2.1 bout file untouched).

Functions (sqw/README.md, "Figure exports"): planAData, planAFlyHMM, planALikelihoods, exportFigA1..A4, planAStyle /
planAFigure / planAStrip. Each export writes the plotted per-fly points and the per-block median + fly-bootstrap CI as
CSV (Origin) and .mat, and draws the figure with the theme pinned light, genotype order CS -> CS DD -> per0 LD ->
per0 DD -> fumin -> sss, light before dark. Outputs: sqw/validation/figA. Tests: test_sqw_figA 7/7 (B = I makes the
HMM marker equal the counted J/EPR to 4e-16; path-mean J vs exact forward-backward posterior current |z| < 0.5;
noiseless Onsager L recovered to 1e-16; exclusion filters, order, theme, Fig A4 clipping and file routing).

**Raw-data checks (same data as stages 1-3; sss = the 54-file 20s-A set for A1-A3, 57 for A4 c/d as 20p/20q):**
- Recipe checks all reproduced: Fig A1 CS dark c_obs +0.147, J +0.167, EPR 0.116; Fig A2 L 0.272 / 0.262 / 0.286 /
  0.352 / 0.243 / 0.145; Fig A3 E_S - E_Q -1.78 (CS dark), +1.00 (fumin dark); Fig A4 dBIC medians -16.6 / -30.4,
  above 0 in 28% / 23%. Every A1/A3 median equals the 20s-A ledger to its printed precision.
- Fig A1 HMM markers (200 FFBS paths per fly, MATLAB stream) against the 20s-A "HMM bracket" medians: 18 comparisons,
  mean z +0.37, max |z| 2.40, none > 3 (SD of the median from 6 independent path draws, so rough). The stored fits
  were matched to data files BY NAME (the fits index the 57-file sss list, 20s-A used 54 files); sss agrees with the
  ledger, so no sign that sqw17 paired sss flies with the wrong fit.
- Fig A4 (a, b): planALikelihoods reproduces the 20n/20o ledger: CS dLL +17.7 [6.4, 30.1], MC2 +36.1, MC1_6 +14.3,
  HMM_6 +31.6, R_mem 2.00 (>= 1 in 0.96), R_drift 0.08 (>= 0.5 in 0.13); fumin +8.7 [1.8, 21.8], +18.7, +13.4, +26.6,
  1.25 (0.57), -0.03 (0.11); dBIC > 10 in 0.20 / 0.21. HMM_6 refitted with 20 MATLAB-stream starts (61 min, 78 flies).
- Fig A4 (c, d) read the per-fly CSVs of the stage-3 validation run (4448209), whose medians were checked there.

**For cClaude (presentation, nothing changed):**
1. fig_beneath.py hard-codes the panel (a) text ("28% (CS), 23% (fmn)", "~3 nats per 1000 min"); the port computes
   it (0.28 / 0.23, 3.1 / 3.0 nats per 1000 min): same numbers.
2. Panel (c): per0 DD is filled (passes the pre-registered rule, a null-based Holm test) while its fly-bootstrap 95%
   CI of the median E crosses 0 (-0.008 to 0.077). Not an error, but a reader comparing the two marks will ask.
3. The 20o ledger's "HMM_6 - MC1_6 (gain left after drift)" is over all flies; nothing in the figure uses it.

## 2026-11-01  Figs C8, C10, C12 and the Fig A4(c) legend — code 2902414 (lClaude, matlab-port)

Merged claude/trusting-allen-7nxrsf (68c9c0b) first. Tests: test_sqw_figA 9/9 (adds C10: cosine amplitude and
mean recovered to 1e-17 on a synthetic 20w-B file; C12: dead flies dropped, no-cycle fly-phases and JSON nulls out of
the depth and J panels, ledger stars as given).

- **Fig A4(c) legend** now carries cClaude's wording (filled = pre-registered test, each fly vs its own simulated
  null; bars = descriptive bootstrap CI of the median, can cross 0 when filled). The x range is widened on the right
  so the four-line note does not sit on the sss rows. Re-rendered from the repo copies of its inputs
  (sqw/validation/figA/figA4_likelihoods_perfly.csv, sqw/validation/perfly_20pq).
- **Fig C8** (analysis_20r/figC8_wake_hazard.m; outputs analysis_20r/figC8): per0 DD (replication bouts, b14ba340,
  n 4017) and CS DD (80402183, n 3309) wake hazards at 12/16/20/24 bins, recomputed with rise_test_20r's settings and
  999 null replicates. All eight RiseNullP equal the ledger (per0 DD 0.001 / 0.009 / 0.059 / 0.086; CS DD 1.000 /
  1.000 / 0.994 / 0.140), and the troughs, peaks and RiseRatios equal the ledger entries. Common log axes in minutes.
- **Fig C10** (sqw/exportFigC10; outputs sqw/validation/figC): from planB_20wB.json; the port refits the 24-h cosine
  and its amplitude equals the stored AJ for all four genotypes (max difference 4e-17). Verdict words copied.
- **Fig C12** (sqw/exportFigC12): from y_kir.json; flies 19 / 25 / 34, depth and J n = 33 for vnc (one no-cycle
  fly-phase each). Medians reproduce the 20y ledger: sleep light 38.5 -> 23.9 (vnc), dark 42.6 -> 34.7 (dfb) / 25.3
  (vnc); E31 raised by 0.37 / 0.58 (dfb light / dark) and 1.08 / 1.05 (vnc). The Holm stars copied in fig_kir.py
  agree with the 20y ledger in every panel. The clipped depth value is found and named from the data (control dark at
  7.7), not hard-coded.

**Defect found in the bout suite, not changed (root rule: ask first):** BoutHazard's optional SurvivalHandle overlay
evaluates S(lefts) with lefts a ROW; FitHyperexponentialMLE's SurvivalHandle returns a COLUMN, so
(sl - sr) ./ sl ./ width broadcasts and H.ModelHazard comes back n x n (11 x 11 on per0 DD, 12 bins). It touches only
the overlay (H.ModelHazard and the model curve drawn by Plot=true); RiseNullP goes through NullSurvivalHandle and the
C8 reproduction is exact. test_BoutHazard test 5 misses it because its handle preserves the input shape. Fix: one
reshape in BoutHazard plus a test with a column-returning handle. figC8 computes the null hazard itself meanwhile.

## 2026-11-01  BoutHazard overlay fixed; Fig A1 observed-only; Fig C7 observed — code cd2f86f (lClaude, matlab-port)

Done on Sheyum's request ("port them"; "fix the error you found in BoutHazard.m") after cClaude's a450c0a.

- **BoutHazard (root bout suite) fixed.** The SurvivalHandle overlay now takes the handle's output as a row, as the
  null sampler in the same file already did. New test_BoutHazard test 11 (a row- and a column-returning handle must
  give the same 1 x nb ModelHazard, equal to the exact rate) fails on the old code with [11 11] and passes now;
  test_BoutHazard 11/11. On the per0 DD bouts the fixed overlay equals figC8's independent calculation to 3e-15 at
  12/16/20/24 bins. RiseNullP never went through this path; the C8 and 20r numbers are unaffected. Root CLAUDE.md
  counts updated (test_BoutHazard 11, 130 bout tests). The analysis-folder copy of BoutHazard.m (Icode_sync_stamp.txt)
  is now behind; it was NOT touched.
- **Fig A1** has no HMM markers by default (main text observed-level only); HmmMarkers=true writes the SI version
  (figA1_direction_flux_SI_hmm). Points and medians unchanged (figA1_points.csv / _summary.csv byte-identical).
- **Fig C7 observed** (exportFigC7; fig_R1_obs.py): routeDurationObserved / routeDurationD ported from sqw_ref and
  checked exactly against golden G11 (47 records from the G2 trace; three D cases, diff 0). The 20k observed R1
  statistics recomputed from the raw records, on 20k's flies (non-degenerate 20b fits: CS 25/25, fumin 53/53,
  per0 LD 27/27, per0 DD 23/25, sss 54/57): n and mean D equal the ledger in all 9 blocks (CS L +0.087, CS D -0.293,
  fumin -0.002 / -0.061, per0 LD -0.228 / -0.121, per0 DD -0.028, sss +0.010 / +0.008); memoryless reference 0.2016
  (recipe 0.202).

**Case (Monte Carlo, nothing changed in the ledger):** per0 LD light is Holm-significant in the ledger (p 0.005,
Holm 0.04) because, as the second-smallest of 9 p-values, its Holm value is 8p: the star needs p < 0.00625. This
port's own 2000-permutation run (seed 20) drew p 0.009, Holm 0.072, i.e. no star. Diagnostic only (not written to
the repo): 20000 permutations, seeds 20 and 21, gave p 0.0041 / 0.0043, Holm 0.033 / 0.034, so the ledger's star is
the right one and the port's draw was high (9 tail permutations where ~4 were expected). The figure therefore takes
its stars and fading from the ledger's Holm by default (SigFrom="ledger", as Fig C12 copies its flags); figC7_blocks
carries both. Suggestion for cClaude: 2000 permutations leave this verdict at the mercy of the draw; if 20k is ever
re-run, more permutations would settle it.

## 2026-11-02  Fig A1 sss mean marker, Fig A2 sss note and panel (c) A_eff (lClaude; FIGURE_RECIPES notes 2026-11-02)
- exportFigA2 panel (c): A_eff = EPR/|J| medians reproduce fig_planA.py's printout exactly: CS 0.57, CS DD 0.51,
  per0 LD 0.43, per0 DD 0.64, fumin 1.11, sss 0.70 nats/cycle (n 48, 25, 51, 24, 88, 76 fly-phases); L unchanged.
  sss title notes '(open = light: reversed)'; top scale e^A (odds forward : backward). Output figA2_Aeff.csv.
- exportFigA1: sss mean (diamond, fly-bootstrap 95% CI, left of the median) in (a) and (b), main and SI variants.
  Values (no Python printout to compare): c_obs light -0.111 [-0.187, -0.042], dark -0.068 [-0.150, +0.001]; count J
  light -0.037 [-0.061, -0.013], dark -0.005 [-0.030, +0.021] cycles/h. Summary CSV gains mean, mean_lo, mean_hi.
- fig_sss_reversal (exportFigSssReversal): against PLANS C2 (fig_sss_reversal.py printout): sss light n 50, chain J mean
  -0.0377 (-0.038), p 0.0045 (0.004); HMM J -0.061 (-0.062) and p 3e-5 (7e-5) are FFBS Monte Carlo. Two deterministic
  numbers differ: J < 0 in 62% (Python 64%) and sss dark p 0.63 (0.59). CASE for cClaude: 13 of 50 sss-light flies
  (5 of 53 dark) have no cycle in their observed chain, so J_chain is round-off (+-1e-15) and its SIGN differs between
  numpy and MATLAB; the "% J < 0" and, through tiny Wilcoxon ranks, p depend on it. Counting |J| < 1e-12 as 0: sss
  light J < 0 in 52% (26/50), p 0.0016; sss dark 51%, p 0.67 (summary columns neg_obs_cycleOnly, p_obs_zeroed).
  The reversal conclusion stands either way; the printed 64% is not a stable number. Figure text mirrors the Python.
- fig1_hmm_emissions (exportFig1Emissions): n 25 / 27 / 23 / 53 / 54 non-degenerate fits; diagonals Q 0.80-0.94,
  W 0.91-0.96, S 1.00 (recipe ranges reproduced: Q 0.916 0.838 0.936 0.831 0.802; W 0.910 0.935 0.912 0.962 0.946).
- fig_dwell_direction_CS_night (exportFigDwellDirection): reproduces the ledger entry exactly (deterministic): 25 CS
  flies; Q into / out of sleep n 649 / 433, per-fly mean 2.68 vs 2.38 min (medians over flies), into longer in 64%;
  W into / out n 418 / 585, 3.47 vs 4.43 min, out longer in 80%.
- C13 profiles fig_planB21_profiles_primary / _supp (exportFigC13Profiles; fig_planB_21.py): against the Python's
  stored fig_planB21_profiles.json, all 12 genotype x setting rows agree: fly counts equal (L300: 25 25 25 27 51 45;
  L120: 25 25 25 27 51 46) and m_start, m_end and median d agree to <= 6e-16 at every bin. Bands are bootstrap.
  per0 DD traces cut to the 20r replication traceLen, as planB_21.traces. Check table fig_planB21_check.csv.
- SI memory figure with HMM panel (fig_R1_main.py): NOT finished at the time of this line; completed in the next entry. routeDurationHidden ported (golden G11 exact;
  test_sqw_figA 11/11) and planARouteHidden written (per-fly hidden R1 D over 200 FFBS paths, 400-simulation null with
  20 paths each, path records, landscape line), smoke-tested on 2 CS flies with 10 sims (26 s). At 400 sims the full
  run is estimated at many hours on 2 workers; not started. Still to do: the full run, exportFigC7SI (hidden panel,
  hidden Holm), and validation against the 20k ledger's hidden R1 lines (D, null mean, p, Holm).
- SI memory figure figS_C7_memory_with_hmm_check (exportFigC7SI, planARouteHidden; fig_R1_main.py), completed:
  hidden-level 20k R1 recomputed from the stored fits (182 non-degenerate flies; 200 FFBS paths per fly; null = 400
  records per fly simulated from its fit, unmasked, each decoded with 20 paths; ~56 min on 2 workers). Against the
  ledger's hidden R1 lines: verdicts identical (only CS dark significant: p 0.0050 = the 400-simulation floor, Holm
  0.0449, as the ledger); mean D within 0.006 in all 9 blocks (CS D -0.280 vs -0.275); n equal in 8 of 9 (sss dark 22
  vs 23: one fly at the "finite in half the paths" threshold, which flips with the path draw); null means and p
  within Monte Carlo error (largest: fumin dark p 0.22 vs 0.35). Stars from the ledger by default (SigFrom), port
  values in figS_C7_hidden_blocks.csv. Panel (c): landscape prediction 0.349 (mean over CS flies of (F_SW+1)/2).
  The per-fly cache (550 MB) stays in the session scratchpad, not the repo.

## 2026-11-02  fig_sss_reversal: no-cycle rule adopted (cClaude's decision, FIGURE_RECIPES 2026-11-02)
- exportFigSssReversal now sets |J_chain| < 1e-12 (20s-A) to exactly 0, plots it at 0 and leaves it out of the
  observed "% J < 0" and of its Wilcoxon test; the panel text prints "(n no cycle)"; the bar is over all flies; HMM J
  untouched. Summary columns now n, n_noCycle, mean_obs, mean_hmm, neg_obs_cycleOnly, neg_hmm, p_obs_cycleOnly, p_hmm
  (the round-off-sign neg_obs / p_obs and p_obs_zeroed are gone).
- Against the corrected Python (RESULTS_LEDGER correction 2026-11-02): sss light 37 of 50 with a cycle, J < 0 70.3%
  (70%), p 0.00161 (0.0016), mean -0.0377 (-0.038) -- match. sss dark 48 of 53, 56.3% (56%) -- match; p 0.675
  against 0.66. CASE for cClaude (not tweaked): sss dark has one exact tie in |J|, two flies at J = -1/24 and +1/24
  (-0.041666666666664 / +0.041666666666664: 11_5_2018_4b_WC_pruned_11 and 3_23_2018_4boardmb_pruned_22), so the
  p depends on tie handling: MATLAB signrank exact with mid-ranks 0.675; normal approximation with tie correction 0.670;
  tie broken either way, exact 0.651 / 0.681. Python's 0.66 is in that range but I cannot tell which rule scipy used
  (no Python or s_observable.json here); its full-precision p and scipy version would settle it. Far from 0.05 either
  way, so the verdict is unaffected. (sss light has no exact tie; closest pair 1.8e-14 relative, p matches.)
- CS for reference: light 24 of 25 with a cycle, 33% negative, p 0.15; dark 25 of 25, 4%, p 3e-7.

## 2026-11-02  Rule of three for one-way edges (convention change, Sheyum) — PORTED FROM THE SPEC, NOT YET VALIDATED
- sqw/ruleOfThreeMatrix.m, sqw/modelEnergeticsR3.m written from the spec relayed by cClaude (sqw_ref section 16 is NOT
  on claude/trusting-allen-7nxrsf or any remote branch as of a33a73f, nor are sqw_r3.py / sqw_out/r3_model.json).
  test_sqw_r3: 8/8 (hand example [180 3 20]/203, min(3, N(b,a)), both-way-unseen pairs kept 0, order independence,
  n_a = actual row steps for an empirical P, finite EPR / A where modelEnergetics gave Inf).
- Interpretations to confirm against the Python: (i) all of a row's one-way edges are filled together from the
  ORIGINAL P (order-free); (ii) the count n_a P(a,b) of the unseen edge (< 1e-10 n_a) is kept, not zeroed, before
  adding min(3, N(b,a)). Corners the spec does not cover, REPORTED not floored (info.stillOneWay, EPR stays Inf):
  N(b,a) = 0 while P(b,a) >= 1e-10 (possible for the marginalised chain fit), and a row with P(a,a) = 1 (n_a = 0/0).
- Not yet wired into observableQuantities / planAFlyHMM / energyLandscape: waiting for the reference to see exactly
  which N each caller uses (observed counts per phase; mean posterior-path counts).
