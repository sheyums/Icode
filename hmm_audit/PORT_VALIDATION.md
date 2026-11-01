# MATLAB port validation (sqw_ref.py -> Icode/sqw), branch matlab-port

Append-only, one entry per stage, with the commit of the code that passed. Port by lClaude.

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
