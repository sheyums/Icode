# PRESPEC 20t — Plan B feasibility checks F1-F6 (before running anything)
Date 2026-10-31. Branch `planB-feasibility` (independent session). Results go to `RESULTS_LEDGER_planB.md` on this branch,
never to the shared ledger. Each check is time-boxed and ends in PASS / FAIL with the number behind it. No check is
re-run with a changed criterion, grid, window or effect size after its result is seen; anything added later is labelled
POST HOC.

## The model (PLANS.md, Plan B), made concrete
Observed 1-min S/Q/W symbols y_t (sqw_ref conventions: Q 0, W 1, S 2, missing -1). Rates per minute
  log k_ij(t) = c_ij + alpha_ij cos(w t) + beta_ij sin(w t) + b_ij H(t),   w = 2 pi / tau,  t in minutes from ZT0,
i.e. a_ij = sqrt(alpha^2 + beta^2), phi_ij = atan2(-beta, alpha) (the linear form makes the fit well behaved).
Per-minute transition probabilities by competing risks (exact for rates constant within the minute, one jump per
minute, the same resolution as the data): K_i = sum_j k_ij, P_ii = exp(-K_i), P_ij = (1 - exp(-K_i)) k_ij / K_i.
Transition into minute t uses rates evaluated at t-1 (cos, sin, H known before the step).
Process S: H in [0,1]; per minute, if y_{t-1} = S: H_t = H_{t-1} exp(-1/tau_d); else (Q, W or missing wake minute)
H_t = 1 - (1 - H_{t-1}) exp(-1/tau_r). Initial value: run the recurrence over the record's first 1440 minutes from
H = 0.5 and start t = 0 at the end value (same rule for truth and fit; the simulator instead discards a 2-day burn-in).
Fit: maximum likelihood, rows independent given H; L-BFGS-B with analytic gradient, started from the constant row.
(tau_r, tau_d) profiled on the FIXED grid {1, 2, 4, 8, 16, 32} h x {1, 2, 4, 8, 16, 32} h (36 points), common to the
three rows of a fly. tau (period) known: 24 h (synthetic). Missing minutes do not contribute transitions.

Per-row tests (per fly): CLOCK in row i = LRT full vs alpha = beta = 0 on both destinations of row i (4 df);
HOMEOSTAT in row i = LRT full vs b = 0 on both destinations of row i (2 df), each re-profiled over the grid.
Nominal chi-square 5% critical values (9.49, 5.99). The chi-square is NOT trusted: (tau_r, tau_d) is unidentified
under b = 0 (Davies problem), so the actual false-positive rate is MEASURED (unmodulated W row; null scenario N0).

## F1 — synthetic recovery of a_ij, b_ij at our data size (synthetic, no data)
Truth (per fly, 25 flies per scenario, seeds fixed: scenario s uses default_rng(20_000 + s)):
- Baseline rates: the G12 golden first-order chain (golden/G12_planA.json, mc_fit, light phase; leaving rates
  Q 0.29, W 0.10, S 0.046 per min, inside the 0.03-0.4 range of the fitted fly chains), converted to rates by
  inverting the competing-risks map; per fly each log rate jittered N(0, 0.25^2). c_ij set so that the record-average
  log rate equals the baseline (the H term is centred on the calibration mean of H).
- Modulated transitions (biologically motivated; the W row is left UNMODULATED as the in-fly null):
  clock on Q->S (peak ZT18) and S->W (peak ZT6); homeostat on Q->S (b > 0) and S->W (b < 0). Other transitions constant.
- Effect sizes, anchored to the observed light/dark contrasts of the fitted chains (log rate ratio dark/light
  0.45-0.6 on Q->S, S->W, W->Q, i.e. a cos amplitude of about 0.35-0.47):
  PRIMARY: a = 0.40; homeostat effect E_b = b * SD_t(H) = 0.20 (comparable to a / sqrt(2) = 0.28).
  SECONDARY (scenario DD-CH-half): a = 0.20, E_b = 0.10.
  b is fixed per scenario from SD(H) in 10 calibration flies (seed 19_999, discarded) so that E_b hits the target.
  Per-fly jitter: amplitudes x U(0.75, 1.25), clock phases + N(0, 1 h), tau_r, tau_d x exp(N(0, 0.2^2)) around
  (6 h, 3 h) (deliberately OFF the fitting grid).
- Scenarios: N0 null (no modulation), 4 d | DD-CH clock + homeostat, 4 d (CS DD analogue; PRIMARY) |
  DD-CH-half (secondary effect size), 4 d | DD-C clock only, 4 d | DD-H homeostat only, 4 d (per0 DD analogue) |
  LD2-CH clock + homeostat, 2 d (shortest LD record) | LD3-step-C DIAGNOSTIC: clock + a light step (+0.4 on Q->S and
  -0.4 on S->W in the dark half, i.e. masking), b = 0, fitted with the cos + H model (is masking absorbed as a spurious
  homeostat?) | N0-2d null, 2 d. Scenario numbers s for the seeds: N0 0, DD-CH 1, DD-CH-half 2, DD-C 3, DD-H 4,
  LD2-CH 5, LD3-step-C 6, N0-2d 7. The full model is fitted in every scenario. Criterion (2) for LD2-CH uses N0-2d.
PASS criteria for F1 (all must hold, in DD-CH AND LD2-CH, primary effect size):
 (1) every truly modulated (row, regulator) term detected (LRT at its nominal 5%) in >= 80% of flies;
 (2) false positives <= 10% for clock and for homeostat in the unmodulated W row, and in N0 for every row;
 (3) cross-talk <= 10%: homeostat detected in Q or S row of DD-C, clock detected in Q or S row of DD-H;
 (4) b within tolerance: per modulated transition, median over flies of E_b_hat / E_b_true in [0.75, 1.25]
     (E_b_hat = b_hat * SD_t(H_hat) at the profiled grid point; scale-free, so it does not depend on tau_r, tau_d being
     recovered) AND sign of b_hat correct in >= 80% of flies.
 Reported, not pass/fail: DD-CH-half detection rates; a_hat / a_true and phase error; profiled (tau_r, tau_d) vs truth;
 destination-level detection; LD3-step-C spurious-homeostat rate.

## F5 — synthetic recovery of the housekeeping / driving split (synthetic, no data)
Per fly, with the rates K(t) taken as the protocol (H along the realised path; see the caveat below), p(t) from the
master equation dp/dt = p K(t), exact per minute (matrix exponential), started at the stationary pi(K(0)); averaged over
the record, nats/h (Esposito & Van den Broeck 2010 PRE 82:011143; Hatano & Sasa 2001 PRL 86:3463):
  housekeeping (adiabatic) sigma_hk = sum_{i != j} p_i k_ij ln[ pi_i(t) k_ij / (pi_j(t) k_ji) ],
  driving (non-adiabatic)  sigma_na = sum_{i != j} p_i k_ij ln[ p_j pi_i(t) / (p_i pi_j(t)) ],
  total = sigma_hk + sigma_na = Schnakenberg EPR of p(t); pi(t) = stationary distribution of the frozen K(t).
Computed from the TRUE (theta, H) and from the FITTED (theta_hat, H_hat) on the same synthetic records of F1.
PASS (DD-CH and LD2-CH): median over flies of sigma_hk_hat / sigma_hk_true in [0.80, 1.25]; median of
sigma_na_hat / sigma_na_true in [0.5, 2.0]; and the noise floor, median fitted sigma_na in N0, is < 25% of the median
true sigma_na in DD-CH. Also reported: sigma_na / total (how big the driving share is at all).
CAVEAT stated now: H is computed from the fly's own sleep, so the homeostat is FEEDBACK, not an external protocol;
the split above conditions on the realised path. Whether that is the right thermodynamic object is a modelling question
F5 cannot answer; F5 only asks whether the defined quantity is recoverable.

## F2 — separability of clock and homeostat in real flies (needs data)
Per fly, H(t) from the fly's own symbols at the reference (tau_r, tau_d) = (6 h, 3 h) PRIMARY, and at every grid point
(range reported). r = multiple correlation of H on {cos wt, sin wt} (= max over clock phase of |corr|; the clock phase is
free in the model, so this is the relevant collinearity). tau = 24 h (LD), per-fly free-running period (chi-square
periodogram 20-28 h, 1-min activity) for CS DD; per0 DD reported with tau = 24 h for reference only.
PASS per genotype if r < 0.8 in more than half of the flies (PLANS "most"); reported fraction and median r. The same r
is reported for the F1 synthetic flies (DD-CH) for comparison.

## F3 — informativeness of the clock-precision bound (needs CS DD data)
Per cycle marker: the circular mean (centroid) time of activity counts within each successive 24-h window
(t = 0 .. 4 d). Intervals between successive markers; per fly its mean interval removed (between-fly period differences
are not jitter). Pooled s^2 = sum of squared deviations / (n_intervals - n_flies). s overestimates intrinsic jitter
(marker noise adds), so B = 2 tau / s^2 is still a valid, weaker lower bound. Fly bootstrap (2000) CI for s.
Rhythmic flies: periodogram peak above the 0.01 Sokolove-Bushell threshold. PASS if >= 15 of 25 CS DD flies are
rhythmic AND the 95% CI of the bound B = 2 tau / s^2 spans less than a factor of 4 (upper/lower). per0 DD (no rhythm)
reported as the negative control (expected: no rhythmic flies, bound 0).

## F4 — deprivation data (needs data; format, n, spread only)
PASS if the files are readable into 1-s or 1-min counts with a known deprivation window, n >= 15 deprived flies, and the
per-fly sleep lost (baseline-day sleep minus deprivation-day sleep in the deprivation window) has a 90th/10th percentile
ratio >= 2 (enough spread for a dose-response). No rates are fitted.

## F6 — real-data pilot, CS DD vs per0 DD (needs data; pre-registered here)
Symbols exactly as sqw_ref.sqw_symbols(theta 0.5, light_dark=False, whole-record reference). Models per fly:
CS DD full (clock at the fly's periodogram tau + homeostat); per0 DD homeostat only, and per0 DD full at tau = 24 h as
the clock NEGATIVE CONTROL. Constant chain = one 3x3 competing-risks chain (= MC1 whole record).
(a) Per fly, LR(full vs constant), p by parametric bootstrap: 19 records simulated from the fly's constant chain, each
    refitted through the whole grid procedure; p = (1 + #{LR* >= LR}) / 20.
(b) Per-row CLOCK and HOMEOSTAT LRTs as in F1, thresholds = the 95th percentiles measured in F1's N0 null at 4 d (the
    sample-size-matched null), not chi-square.
(c) 20q follow-up: S1_mod = B_obs / median(B*) with B* from 50 records simulated from the fly's FITTED MODULATED model
    (sqw_ref.tur_estimate, tau = 3 h, same estimator and mask as 20q); compare with 20q's S1 (CS DD 0.75, per0 DD 0.92).
PASS if: modulated model beats constant (p < 0.05) in >= 50% of CS DD flies; CS DD clock term detected (row test) in
>= 50% of flies in at least one row; per0 DD clock negative control detected in <= 15% of flies in every row; and
the 20q extra variance "vanishes": CS DD median S1_mod in [0.9, 1.1] with Wilcoxon of log S1_mod vs 0 p > 0.05.
F6 is run only after F1 passes (otherwise its row tests cannot be read).

## Step-4 decision (PLANS.md) applied as written
Merge only if F1, F2 and F6 all PASS and a figure-level Plan B result exists without reopening Plan A claims;
otherwise paper 1 = Plan A + outlook, Plan B = paper 2. A check that cannot run for lack of data is NOT a pass.

Narrowest block for every real-data number: one fly's whole DD record (F2 LD: one fly's whole LD record).
Sample-size-matched null: N0 at the same record length (4 d DD; 2 d LD for LD2-CH, which uses an N0 run at 2 d as well,
seed 20_007, added to the scenario list now).
