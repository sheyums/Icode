# PRESPEC 20 — 3-symbol (S/Q/W) per-minute HMM: routing through H1 with distinct emissions (before running)
Date 2026-10-23. Requested by the user. Motivation: in the scored 3-state HMM, H1 emits exactly like H2 or H3 in
most windows, so routing (Xn, Yn) and energetics are not identified (PRESPEC 1, 2). Giving quiet wake its own
observable symbol may give H1 a distinct emission row, which is what identification needs.
Caveat fixed now: the Q/W threshold reuses PRESPEC 19's theta and reference, chosen after the +-5 min inspection
(PRESPEC 18). Every p-value here is DESCRIPTIVE. The threshold-free version is a separate, later plan (TODO.md).

## 1. Observable sequences (1-min steps, fixed clock bins)
- Data: count files, all 5 groups (CS, fumin, per0LD, per0DD, sss). Minute m = seconds [60m, 60m+60) from ZT0.
- Sleep = inactive run >= 300 s with bridging (b) (`v18.sleep_runs(x, bridge=True)`).
- Minute symbol:
  - **S** if >= 30 of its 60 s lie inside a sleep bout;
  - otherwise a wake minute with activity a_m = fraction of its seconds with count > 0:
    **Q** if a_m < theta * r_f, else **W**.
- r_f = the fly's reference activity, as in PRESPEC 19: median over that fly's wake bouts >= 10 min in the block of
  the active fraction in the bout's middle 5 min. Blocks: light / dark (ZT 0-12 / 12-24); per0DD whole record.
- theta = 0.5 PRIMARY; 0.33 and 0.75 sensitivity, all reported.
- Symbol codes (for the .mat export and the B matrix columns): 1 = Q, 2 = W, 3 = S.

## 2. Windows
- Fixed 4 h clock windows: ZT 0-4, 4-8, ..., 20-24 on every day (240 steps each). **Never extended to finish an
  event** (length bias, PRESPEC 8b). Incomplete trailing windows are dropped.
- A window is FIT if it contains >= 1 S minute and >= 1 wake minute; otherwise it is reported as unfit.

## 3. Model fits per fly-window: K = 2, 3, 4 with BIC
- Discrete-emission HMM (3 symbols), EM (hmmlearn CategoricalHMM), 20 random starts per K, seed 20 + window
  index; best log-likelihood kept. Initial distribution free.
- Free parameters k = (K-1) + K(K-1) + 2K: k = 7, 14, 23 for K = 2, 3, 4. n = 240 (number of 1-min steps).
- BIC = -2 logL + k ln n. Record logL and BIC for every K.
- **Flag rule (fixed now):** a window is "strongly not 3-state" if K = 2 or K = 4 beats K = 3 by > 10 BIC units.
  Flagged windows are REPORTED, not dropped.
- Reported: per genotype x block, the fraction of windows best at K = 2/3/4, the distribution of
  BIC(K) - BIC(3), and the modal K. If the modal K over all windows is not 3, a secondary fit of every window at
  the modal K is run and labelled SECONDARY; routing in the paper's sense is defined only at K = 3.
- Known bias, stated now: BIC charges for parameters the data cannot identify, so it leans toward fewer states.

## 4. Routing at K = 3 (primary)
- State labels from the fitted emissions: H3 = highest P(S); of the other two, H2 = highest P(W); remaining = H1.
- **Identifiability check per window:** H1 is degenerate if its emission row lies within total variation 0.05 of
  H2's or H3's row. Degenerate windows are reported and excluded from Xn/Yn summaries (their routing is not
  determined, PRESPEC 1).
- Predicted routing per window from P, same formulas as `common.py`: Xn_pred, Yn_pred.
- **Route counts from whole-record decoding**, not window extension: one time-inhomogeneous HMM per fly in which
  each minute uses its own window's (label-aligned) parameters. FFBS, 200 posterior paths. On each path,
  collapse runs to a sequence of distinct states and count N312, N32 (sleep -> wake via H1 / direct) and N213,
  N23 (wake -> sleep via H1 / direct). Each route is attributed to the window in which it COMPLETES (arrival in
  the destination state). Viterbi counts computed the same way, for comparison only.
- Xn = (N312 - N32)/(N312 + N32), Yn = (N213 - N23)/(N213 + N23) from posterior-mean counts, per fly per block
  (counts summed over the block's non-degenerate windows). Flies need >= 3 routes in each denominator.
- Model-based signed current c_model = (Yn - Xn)/2 (A <-> H2, S <-> H3, Q <-> H1; PRESPEC 19 notes).

## 5. Uncertainty
- Parametric bootstrap on a pre-chosen subset: 6 non-degenerate windows per genotype x block, drawn with seed 20
  before any results are seen. 50 replicates each: simulate from the fitted K = 3 model, refit (20 starts),
  recompute Xn_pred and Yn_pred. Reported: SD and bias. Routing in a window is called "estimable" if SD < 0.2.
- Group summaries: mean Xn, Yn, c_model per genotype x block with a 95% fly-bootstrap interval (2000, seed 20).

## 6. Energetics per K = 3 window (non-degenerate only)
From the minute-step P and its stationary pi:
- EPR = sum_{i<j} (pi_i P_ij - pi_j P_ji) ln(pi_i P_ij / (pi_j P_ji)) per minute, converted to per hour;
- cycle flux J = pi_1 P_12 - pi_2 P_21;
- affinity = ln(P_12 P_23 P_31 / (P_13 P_32 P_21)).
A zero transition probability makes EPR and affinity infinite. Such windows are counted and reported, never
floored. Group medians with fly-bootstrap intervals.

## 7. Synthetic validation (run FIRST; the real data are not fitted until it is logged)
- **V1, positive control.** Truth: a 3-state model with distinct H1 emissions and known Xn = -0.3, Yn = +0.4,
  dwell times like the data's; 50 windows x 240 steps.
  Pass: K = 3 chosen by BIC in >= 70% of windows, and |mean bias| < 0.1 in Xn_pred and Yn_pred.
- **V2, negative control.** Truth: a 2-state model (S/W with Q emitted by wake at the data's rate), no H1.
  Pass: K = 2 chosen in >= 80% of windows, and the K = 3 fit's H1 is flagged degenerate, or its routing has
  bootstrap SD >= 0.2, in >= 80% of windows.
- Failure of either is reported in the ledger, and the design is reviewed before touching the fly data.

## 8. Comparison with the model-free current (PRESPEC 19)
- Per fly, c_model (section 4) vs model-free c at the same theta: Spearman correlation per genotype x block.
- Prediction stated now: the group mean c_model has the same sign as the model-free c in each genotype x block
  where the model-free interval excludes 0 (CS dark, fumin light/dark, per0DD, sss light).

## 9. Outputs
- `sqw.py` (build + fit + routing), `sqw_synth.py` (section 7), `sqw.json`, `sqw.log`, ledger entries.
- `.mat` export per genotype in the uploaded layout: `R{fly, window}` = 1 x 240 symbols (1 = Q, 2 = W, 3 = S);
  `C{3*fly-2, window}` = model index (3), `C{3*fly-1, window}` = P (3 x 3), `C{3*fly, window}` = B (3 x 3).
  Files go to the scratchpad, NOT the repo (derived from the user's raw data).
- Figure: three panels on identical Xn-Yn axes [-1, 1]: (a) observed per fly from counts,
  Xn_obs = 2 P(Q|waking) - 1, Yn_obs = 2 P(Q|onset) - 1; (b) S/Q/W HMM, FFBS; (c) S/Q/W HMM, Viterbi.
  Plus a BIC panel (fraction of windows best at K = 2/3/4 per genotype).

## Amendment A (2026-10-23, before any fitting; prompted by the user's question on BIC's n)
With n = 240 the BIC penalty is the same in every window whatever its number of bouts, but transition parameters
are informed by state CHANGES, not steps. So:
- Primary rule unchanged (n = 240).
- Sensitivity rule, pre-specified: n_eff = number of observed symbol changes in the window (minimum 2),
  BIC_eff = -2 logL + k ln n_eff. Report how often the chosen K differs between the two rules.
- The number of observed symbol changes is reported next to every window's BIC values; Delta-BIC is not
  compared across windows as strength of evidence.
- Parametric-bootstrap LRT on flagged windows only (Delta-BIC > 10 rule): up to 20 flagged windows per
  genotype x block, chosen with seed 20; 100 replicates simulated from the K = 3 fit; each refit at K = 3 and at
  the competing K with 10 starts; p = (1 + #{LR_rep >= LR_obs}) / 101. Descriptive.
- The same n_eff variant is computed in the synthetic validation (section 7), reported but not a pass criterion.
