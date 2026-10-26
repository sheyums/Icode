# Figure recipes: data extraction and plotting, for a from-scratch MATLAB implementation

Audience: the Claude instance or person writing MATLAB functions that rebuild these figures. Read
`PORTING_GUIDE.md` first; it defines the conventions used here (symbol codes, state labels, phases, signs,
Jensen, rule of three). Function names in `code` are from `sqw_ref.py`, the reference implementation. Each recipe
ends with numbers from `RESULTS_LEDGER.md` that a correct MATLAB version must reproduce.

Common inputs:
- 1-s count files, one per fly (`counts/<group>/**/*.txt`), groups CS, fumin, per0LD, per0DD (constant
  darkness: no light phase), sss.
- Fitted 12 h models, one per fly (PRESPEC 20b): pi (3), P (3 x 3 x 2: light, dark; one matrix for per0DD),
  B (3 x 3), labelled H1, H2, H3; flies whose H1 is degenerate are excluded.
- theta = 0.5 unless a recipe says otherwise.

Common plotting style (all figures):
- Fonts: STIX (Times-like); 8 pt; axes line 0.6 pt; ticks inward.
- Widths: 7.0 in (double column) or 3.4 in (single).
- Genotype colours, fixed order (validated for colour-vision deficiency): CS #2a78d6, fumin #eb6834,
  per0LD #1baf7a, per0DD #eda100, sss #e87ba4.
- Light phase: open markers; dark or DD: filled markers.
- Never encode identity by colour alone: direct labels or a legend are always present.

---------------------------------------------------------------------------------------------------------------

## A. Observed runs and dwell times (Fig. 9, data only; Fig. 9s)

### A1. Build the minute sequence for one fly
1. Read counts `x`.
2. Score sleep bouts: inactivity >= 300 s, with lab bridging (b) (`sleep_bouts(x, bridge=True)`; see
   PORTING_GUIDE section 4).
3. Make S/Q/W minute symbols `y` (`sqw_symbols(x, theta, light_dark)`):
   - S if >= 30 s of the minute lie in a sleep bout;
   - otherwise Q if the minute's active fraction < theta * r_f of its phase, else W.
   - r_f = median, over wake bouts >= 10 min whose waking or onset falls in that phase, of the active fraction in
     the bout's central 5 min. At least 3 bouts are needed; otherwise that phase's wake minutes are missing (-1).
4. Truncate to whole days.
5. Phase per minute: light if mod(minute, 1440) < 720 (ZT 0-12); per0DD has a single phase.

### A2. Extract runs with censoring (`runs_censored(seq, ph, values, nph)`)
- Run types: S = {S}; wake = {Q, W}; Q = {Q}; W = {W}.
- A run is a maximal block of minutes whose symbol is in the set, started at its first minute. It is assigned to
  the phase of that first minute.
- Duration = number of minutes.
- The run is CENSORED (event = 0) if it is cut by:
  - the end of its phase (minute index where the phase changes);
  - a missing minute;
  - the end of the record.
- Otherwise it ENDED (event = 1): the next minute is a valid symbol outside the set.
- When a run is cut at a phase boundary, its continuation in the next phase is NOT counted as a new run (its start
  is unobserved). Skip to the end of that block.
- Exact test: `golden/G8_dwell.json` (runs with and without phase censoring, KM, weighted KM, hazard, geometric).

### A3. Pool and compute survival
1. For each genotype x phase, concatenate (duration, event) over all flies.
2. Kaplan-Meier survival S(t) = P(run > t) (`km_survival`): at each distinct duration t, multiply by
   (1 - events at t / runs with duration >= t).
3. Per-fly median: KM per fly, first t with S(t) <= 0.5.

### A4. Plot Fig. 9 (data only)
- 2 x 4 panels. Rows: light; dark / DD (per0DD in the dark row). Columns: sleep (S), wake (Q or W), quiet wake
  (Q), active wake (W).
- Step plot of S(t), post-step. x log 1-720 min; y log 1e-3 to 1.05. One line per genotype in its colour.
- Caption must state:
  - 1-min grid and theta = 0.5;
  - censoring at phase boundaries;
  - the 5-minute floor of S runs comes from the sleep criterion.
- Check (ledger, PRESPEC 20h): per-fly KM median S run, CS light 14 min, CS dark 27, fumin light 7,
  per0LD light 22; run counts CS light S 1777, wake 1822.

### A5. Fig. 9s (hazards; supplementary)
- Life-table hazard per minute (`life_table_hazard`):
  - bins [lo, hi) with edges 1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144, 233, 377, 721 min (fixed in advance);
  - a run of duration d is at risk at minutes 1..d, so its minutes at risk in a bin = max(0, min(d, hi-1) - lo + 1);
  - events = ended runs with lo <= d < hi;
  - h = events / minutes at risk;
  - do not show bins with < 10 runs at risk.
- Plot at the bin's geometric midpoint, x log. Grid: rows = genotype x phase (9), columns = S, wake, Q, W.
- The version in the repo overlays the fitted-model 95% band (grey); a data-only version omits it.

---------------------------------------------------------------------------------------------------------------

## B. Hidden dwell times (Fig. 10, Fig. 10s)

1. For each fly, sample 200 hidden paths from the posterior (FFBS) under its fitted model, over the whole record
   (`ffbs`; PORTING_GUIDE section 4).
2. On every path, extract runs of H1, H2, H3 with the same censoring rules as A2.
3. Pool runs over paths and flies, each run weighted 1/200, and compute a weighted Kaplan-Meier: weighted events
   over weighted runs at risk.
4. Geometric reference:
   - each fly's dwell survival P(dwell > t) = P_ii^t, from its phase matrix;
   - average over flies, weighted by each fly's number of decoded runs of that state. Equal weighting mismatches
     the pooled KM; this was found and corrected.
5. Mean dwell 1/(1 - P_ii) per fly and phase; plot the median over flies.
6. Fig. 10 layout:
   - 2 rows (light; dark / DD) x 3 columns (H1, H2, H3): decoded KM solid, geometric dashed, log-log axes as in A4;
   - fourth column spanning both rows: median mean dwell per state; x = H1, H2, H3; y log; one marker per
     genotype x phase, open for light.
7. Fig. 10s: log-binned density (edges as in A5) of the weighted decoded dwells, per state and phase.
8. Caption must state that decoded hidden dwells are shaped by the fitted model and are near-geometric by
   construction; they show time scales, not evidence for the model.
9. Check (ledger): H1 mean dwell 2.9-4.0 min in every block; H3 CS dark 35.9 min, fumin dark 8.2 min.

---------------------------------------------------------------------------------------------------------------

## C. Round trips around sleep bouts (Fig. 8) and the per-fly table

### C1. Observed (`round_trips_observed(x, theta, light_dark)`)
1. Use sleep bouts whose preceding AND following wake bouts are each >= 600 s long.
2. Entry is Q if the active fraction in seconds [onset-300, onset-60) < theta * r_f(phase of onset); exit is Q if
   the active fraction in [end+60, end+300) < theta * r_f.
3. Types:
   - forward = (entry Q, exit A);
   - reverse = (A, Q);
   - direct = (A, A);
   - quiet_both = (Q, Q).
4. Count per phase of the onset.

### C2. Hidden (`round_trips_hidden(z, ph, nph, min_sleep, min_gap)`)
1. On the run-compressed path, every H3 run with neighbours on both sides is a round trip. Entry = previous state,
   exit = next state, with H1 standing for Q and H2 for A.
2. Matched definition (used in Fig. 8):
   - H3 run >= 5 min;
   - the previous H3 run ends >= 10 min before it starts, and the next H3 run starts >= 10 min after it ends (both
     must exist).
3. Phase = phase of the H3 run's first minute.

### C3. Per fly (`round_trip_stats(counts, hours)`)
- Fractions of the 4 types, rates per hour, P(quiet entry) = forward + quiet_both, P(quiet exit) =
  reverse + quiet_both.
- Dependence = fraction - product of marginals.
- Hidden level: compute on each path, then average over paths, skipping paths with no round trip.
- Flies need >= 3 round trips in a phase.

### C4. Fig. 8
- 2 panels side by side, (a) observed, (b) hidden matched.
- Rows: CS L, CS D, per0LD L, per0LD D, per0DD, fumin L, fumin D, sss L, sss D.
- Horizontal stacked bars of mean fraction over flies, in the order forward (#2a78d6), reverse (#eb6834), direct
  (#b9b8b3), quiet both (#e4e3df). White 0.8 pt gaps between segments. Percent printed inside segments >= 8%.
- Legend below: "forward Q->S->A (H1->H3->H2)", "reverse A->S->Q (H2->H3->H1)", "direct A->S->A",
  "quiet both Q->S->Q".
- Per-fly table: `roundtrip_table.py` writes `roundtrip_per_fly_theta0.5.csv` (columns listed in
  PORTING_GUIDE section 8).
- Check: CS dark observed forward 0.24, reverse 0.11, direct 0.56; hidden matched forward 0.38, reverse 0.11.

---------------------------------------------------------------------------------------------------------------

## D. Cycle current and entropy production (Figs. 6, 6b)

### D1. Hidden level, per fly and phase
1. 200 FFBS paths; pair counts N_ij on each run-compressed path, by phase of the arrival minute (`pair_counts`).
2. J_cycle = (N32 - N23) / hours of the phase. Forward-positive; linear, so averaging paths or counts is the same.
3. EPR = Schnakenberg sum over pairs with the rule-of-three bound (`epr_count_bound`, identical to the user's
   `epr_count_based.m`, rule 'bound'), on EACH path, then averaged over paths.

### D2. Observable level
- The same counts on the run-compressed S/Q/W symbols (`observable_pair_counts`; Q, W, S in the roles of H1, H2,
  H3; changes touching a missing minute are skipped).

### D3. Group points
- Median over flies (non-degenerate) of J and of EPR.
- 95% bootstrap interval of the median: resample flies, 2000 times, seed fixed.

### D4. Fig. 6 (single panel, 3.4 x 2.9 in)
- x = J_cycle (cycles/h, + = S->W->Reg->S); y = EPR lower bound (nats/h).
- Error bars both ways; marker open = light, filled = dark/DD.
- Dashed grey lines through the origin at A = EPR/J = 0.5, 1, 2 nats/cycle.
- Grey band 0.255-0.327 nats/h: reversible-null 95th percentile, single fly (PRESPEC 20d).
- Grey band 0-0.025 nats/h: the hidden estimator's additive offset (PRESPEC 20f).
- Weak-circulation blocks drawn at alpha 0.35: CS light, per0LD light, sss light, sss dark.
- Labels with leader lines near the origin.

### D5. Fig. 6b
- 2 panels, (a) hidden, (b) observable, same axes x -0.07..0.32, y 0..0.42.
- Same faded blocks; the offset band only in (a).

### D6. Check
- CS dark hidden J median +0.201, EPR 0.301.
- Observable J +0.167, EPR 0.116.
- per0DD hidden +0.233 / 0.208.

---------------------------------------------------------------------------------------------------------------

## E. Three-way Xn-Yn comparison (Fig. 7)

### E1. Columns
- (a) Observable, model-free, whole record per fly (`cycle_current_model_free` counts):
  - Xn = 2 nSQ/(nSQ + nSA) - 1;
  - Yn = 2 nQS/(nQS + nAS) - 1.
- (b) S/Q/W HMM:
  - route counts on each FFBS path summed over both phases (`route_counts`);
  - Xn, Yn per path (>= 3 routes in a denominator);
  - mean over paths, kept if defined on >= 50% of paths.
- (c) Legacy binary model:
  - Viterbi route counts per accepted 4 h window of the scored .mat files;
  - summed over windows per fly, then Xn = (N312 - N32)/(N312 + N32), Yn = (N213 - N23)/(N213 + N23).

### E2. Layout
- 4 x 3 panels; rows CS, per0 (LD green, DD yellow; legacy per0 dark grey #52514e), fumin, sss.
- Equal axes -1.05..1.05, aspect 1; dashed diagonal Yn = Xn; light grey zero lines.
- Text in each panel: n flies and % above the diagonal.

### E3. Check
- % above the diagonal: CS 60 / 88 / 88; fumin 70 / 91 / 100; sss 37 / 32 / 100; legacy per0 44.

---------------------------------------------------------------------------------------------------------------

## F. Things that went wrong once (do not repeat)
- **Weighting the geometric reference per fly while the decoded KM pools runs** (Fig. 10 first draft).
- **Treating the count EPR as a lower bound on the true rate in weak blocks.** It is biased upward there by
  finite counts (PRESPEC 20f).
- **Substituting a flat 3 for a zero count.** It must be b* = min(3, k), in both the prefactor and the log.
- **Comparing hidden round trips (every H3 run) with observed ones (bouts with >= 10 min wake on both sides).**
  Use the matched definition.
- **Reporting nonlinear statistics computed on path-averaged counts** (Jensen). Compute on each path, then
  average.

---------------------------------------------------------------------------------------------------------------

## G. Per-fly EPR vs J_cycle, linear vs quadratic (Fig. 11; PRESPEC 20j)

### G1. Points
- One point per fly x phase:
  - J = count J_cycle averaged over 200 posterior paths;
  - EPR = rule-of-three bound per path, averaged (as D1).
- Exclude:
  - degenerate flies;
  - fly-phases with no cycle (|J_model| < 1e-12);
  - fly-phases without a Q/W reference (PRESPEC 20i).
- Light: open marker. Dark / DD: filled.

### G2. Fits per genotype (both phases pooled), unweighted least squares with a free intercept (`epr_j_models`)
- M1: EPR = b + a|J|.
- M2: EPR = b + cJ^2. c = 1/L in the user's Onsager notation.
- Compare by residual sum of squares (same number of parameters).
- Exponent: least-squares slope of log(EPR - 0.02) on log|J|, using points with |J| >= 0.05 and EPR >= 0.05.
- 95% CIs: bootstrap over flies, 2000 resamples, keeping both of a fly's phases together.

### G3. Layout
- 2 x 3 panels.
  - (a) All genotypes overlaid.
  - (b)-(f) One genotype each, with the M1 fit (solid) and the M2 fit (dashed) over the data range; legend gives a
    and c.
- Title of each panel: n points and the exponent.

### G4. Check
- CS: c = 4.48, exponent 1.76, quadratic better.
- fumin: a = 2.23, linear better.
- sss: c = 6.44.

---------------------------------------------------------------------------------------------------------------

## H. Onsager conductance L per genotype (Fig. 12; post hoc, user request)

### H1. Points
- The same per-fly x phase points as G1: J = count J_cycle, path mean; EPR = rule-of-three bound, path mean.
- Same exclusions as G1.

### H2. Fit (the user's Methods Eq. 25 scheme; `onsager_fit`)
1. Per genotype, regress EPR on J^2 by ordinary least squares THROUGH THE ORIGIN:
   - slope s = sum(J^2 * EPR) / sum(J^4);
   - L = 1/s.
2. Standard error of s: sqrt( RSS / (n - 1) / sum(J^4) ), with df = n - 1.
3. L's standard error: SE_L = L^2 * SE_s.
4. 95% interval for s: s +- t(0.975, n-1) * SE_s. Invert it for L; the bounds swap.
5. R^2 through the origin = 1 - RSS / sum(EPR^2).
6. Do the fit three ways: pooled over phases (main), light only, dark only. per0DD has one phase; plot it with the
   dark points.

### H3. Intervals that respect the fly as the unit (`onsager_bootstrap`, `bootstrap_difference_p`)
- Fly bootstrap: resample flies with replacement, 2000 times, keeping both phases of a fly together; refit L each
  time; 2.5 and 97.5 percentiles.
- Genotype vs CS: bootstrap the difference in L (independent resamples of each genotype, 4000 each); two-sided
  p = 2 * min(P(diff <= 0), P(diff >= 0)).
- The Wald z-test on s, with Benjamini-Hochberg over the 4 contrasts, can be reported alongside. It treats a
  fly's two phases as independent and therefore overstates significance (per0DD vs CS: Wald p ~ 1e-13,
  bootstrap p = 0.065).

### H4. Sensitivity
- Refit with EPR - 0.02 (the calibrated estimator offset, PRESPEC 20f; negative values set to 0).
- L changes by 1-6%; report it, but don't plot it.

### H5. Layout (7.0 x 2.8 in, width ratio 1 : 1.4)
- (a) Pooled phases.
  - x = genotype in the order CS, per0LD, per0DD, fumin, sss.
  - y = L, axis from 0.
  - Filled marker in genotype colour, with fit 95% CI error bars (capsize 2.5).
  - The fly-bootstrap interval as a thin dotted vertical line just to the right of each marker.
- (b) By phase.
  - Light (open marker) at x - 0.12, dark (filled) at x + 0.12; per0DD at x.
  - Error bars = fit 95% CI.
- y label: "Onsager conductance L (cycles^2 nats^-1 h^-1)".
- Footnote: "EPR = J^2/L fitted through the origin per genotype (per-fly points, hidden states). Bars: 95% CI from
  the fit; dotted: fly bootstrap."

### H6. Check (ledger, Onsager L entry)
| genotype | L | fly-bootstrap CI |
|---|---|---|
| CS | 0.215 | 0.171-0.229 |
| per0LD | 0.245 | 0.226-0.271 |
| per0DD | 0.330 | 0.203-0.361 |
| fumin | 0.237 | 0.171-0.274 |
| sss | 0.144 | 0.119-0.178 |

- By phase: CS light 0.297, CS dark 0.213.
- Bootstrap vs CS: sss p 0.008, per0LD p 0.011.
- Code: `onsager_L.py`.
