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
