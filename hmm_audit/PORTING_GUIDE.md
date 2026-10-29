# Porting guide: S/Q/W hidden-Markov analysis, Python -> MATLAB

Audience: the Claude instance (or person) translating this analysis into MATLAB functions for repeated use.
Read this file first, then `sqw_ref.py`. Step-by-step data extraction and plotting for every figure: `FIGURE_RECIPES.md`. The paper's Methods are a separate document (`METHODS_PAPER.md`).

## 0. Sources of truth, in order

1. **`PRESPEC_19.md`, `PRESPEC_20.md` (+ amendment A), `PRESPEC_20b.md` ... `20e.md`.** These define every
   quantity. If code and plan disagree, the plan wins, and the disagreement is a bug to report.
2. **`sqw_ref.py`.** One function per step. It is checked by `verify_ref.py` (log: `verify_ref.log`, ALL PASS)
   to reproduce exactly what the analysis scripts produced: S/Q/W symbols for 187/187 flies; EM
   log-likelihoods to 1e-6; FFBS and Viterbi route and pair counts; sign identities.
3. **`golden/G1..G6.json`.** Fixed inputs with expected outputs, made by `make_golden.py`. Port each function
   until its golden test passes.
4. **`RESULTS_LEDGER.md`.** Every result, with the commit that produced it. A MATLAB port that reproduces the
   fly-level numbers there (given the same fitted parameters) is correct.

Do not port the analysis scripts (`sqw2.py`, `sqw3.py`, ...). They contain paths, multiprocessing and history.
Port `sqw_ref.py`.

## 1. Pipeline

```
1-s counts (txt)                       read_counts
  -> sleep bouts (>= 300 s, bridged)   sleep_bouts
  -> wake-bout events, r_f per phase   wake_bout_events, wake_reference
  -> 1-min S/Q/W symbols               sqw_symbols            (theta = 0.5; 0.33, 0.75 sensitivity)
  -> phase index per minute            phases                 ('12h', '6h', 'none' for DD)
  -> fit K = 2, 3, 4                   em_fit (uses em_step logic), n_params, bic, n_eff
  -> label K = 3 states                label_states, is_degenerate
  -> predicted routing per phase       xn_pred, yn_pred, c_pred, model_energetics
  -> posterior paths                   ffbs (200 paths), viterbi (comparison only)
  -> counts per path                   route_events, route_counts, pair_counts
  -> per-path statistics, averaged     posterior_summaries: Xn, Yn, c, J_cycle, EPR bound
  -> observable checks                 c_obs_minutes, cycle_current_model_free, observable_pair_counts
  -> round trips per sleep bout        round_trips_observed, round_trips_hidden, round_trip_stats
  -> simulation (validation, nulls)    simulate
```

## 2. Function map

Suggested MATLAB names follow the Icode style: CamelCase, an `arguments` block for validation, and
name-value options. Shapes are given for MATLAB with 1-based indices.

| Python (`sqw_ref.py`) | MATLAB (suggested) | In | Out | Golden |
|---|---|---|---|---|
| `read_counts(path)` | `x = readCounts(path)` | text file, 1 integer per line | `x` (T_s x 1) int | -- |
| `sleep_bouts(x, bridge, min_len)` | `bouts = sleepBouts(x, Bridge=true, MinLen=300)` | counts | (nB x 2) `[start, len]` | G1 |
| `wake_bout_events(x, bouts)` | `ev = wakeBoutEvents(x, bouts)` | counts, bouts | struct array (ref, wakePhase, wakeAct, onsetPhase, onsetAct) | G2 (`n_wake_events`) |
| `wake_reference(ev, phase)` | `r = wakeReference(ev, Phase=[])` | events | scalar or NaN | G2 (`refs`) |
| `sqw_symbols(x, theta, light_dark)` | `[y, refs] = sqwSymbols(x, Theta=0.5, LightDark=true)` | counts | `y` (T x 1) in {1=Q, 2=W, 3=S, 0=missing} | G2 |
| `phases(T, scheme)` | `ph = minutePhases(T, Scheme="12h")` | T | `ph` (T x 1) in 1..nP | G3 (`ph`) |
| `emission_lik`, `forward`, `loglik` | `[A, c] = hmmForward(pi, P, B, y, ph)` | params | scaled alpha, scale factors | G3 (`loglik`) |
| `em_step` | `[pi, P, B, ll] = hmmEmStep(pi, P, B, y, ph)` | params | updated params | G3 (`em_step`) |
| `em_fit` | `fit = hmmFitPhased(y, ph, K, Starts=20, MaxIter=1000, Tol=1e-7)` | symbols | struct ll, pi, P (K x K x nP), B, iter | statistical (sec. 5) |
| `n_params`, `bic`, `n_eff` | `k = hmmNumParams(K, nP)`, `b = hmmBic(ll, K, nP, n)` | -- | -- | G3 |
| `label_states`, `is_degenerate` | `[P, B, pi] = labelStates(P, B, pi)`, `tf = isDegenerateH1(B, 0.05)` | K = 3 fit | H1, H2, H3 order | G4 |
| `xn_pred`, `yn_pred`, `c_pred` | `[Xn, Yn, c] = routingPredicted(P)` | 3 x 3 | scalars | G4 |
| `stationary`, `model_energetics` | `[J, EPR, A, piS] = modelEnergetics(P)` | 3 x 3 | per hour | G4 |
| `ffbs` | `Z = hmmFFBS(pi, P, B, y, ph, NPaths=200)` | params, symbols | (NPaths x T) states | G6 (statistical) |
| `viterbi` | `z = hmmViterbiPhased(pi, P, B, y, ph)` | params | (T x 1) | G3 (`viterbi`) |
| `route_events`, `route_counts` | `[names, t] = routeEvents(z)`, `C = routeCounts(z, ph, nP)` | path | events; (nP x 4) | G5 |
| `pair_counts` | `N = pairCounts(z, ph, nP)` | path | (3 x 3 x nP) | G5 |
| `xy_counts` | `[Xn, Yn] = routingFromCounts(C, MinRoutes=3)` | counts | NaN if a denominator < 3 | G5 |
| `j_cycle_count` | `J = jCycleCount(N, hours)` | 3 x 3 | per hour | G5 |
| `epr_count_bound` | use the user's `epr_count_based.m` (rule `'bound'`) | 3 x 3, hours | nats/h | G5 |
| `posterior_summaries` | `S = posteriorSummaries(Z, ph, nP, hours)` | paths | per-path means | G6 |
| `c_obs_minutes` | `c = cObsMinutes(y, ph, nP)` | symbols | per phase | G5 |
| `cycle_current_model_free` | `R = cycleCurrentModelFree(ev, Phase, Theta, hours)` | events | c, J, counts | -- |
| `simulate` | `[z, y] = hmmSimulatePhased(pi, P, B, ph)` | params | path, symbols | statistical |
| `observable_pair_counts` | `N = observablePairCounts(y, ph, nP)` | symbols | (3 x 3 x nP) symbol-change counts; feed to `jCycleCount` / `epr_count_based.m` for observable J and EPR | G7 |
| `round_trips_observed` | `C = roundTripsObserved(x, Theta=0.5, LightDark=true)` | 1-s counts | (nP x 4) counts: forward, reverse, direct, quiet_both | G7 |
| `round_trips_hidden` | `C = roundTripsHidden(z, ph, nP, MinSleep=0, MinGap=0)` | one hidden path | (nP x 4); matched definition: `MinSleep=5, MinGap=10` | G7 |
| `round_trip_stats` | `S = roundTripStats(C, hours)` | (4,) or (paths x 4) | fractions, rates/h, P(Q entry), P(Q exit), dependence (fraction - product of marginals); per path then averaged | G7 |
| `onsager_fit` | `F = onsagerFit(J, EPR)` | per-fly J, EPR | slope 1/L, SE, L, SE_L, 95% CI of L (inverted t interval), R^2 through origin; Eq. 25 | G10 |
| `onsager_bootstrap` | `Lb = onsagerBootstrap(points, NBoot=2000)` | (fly, J, EPR) rows | bootstrap L, flies resampled with all their points | G10 (statistical) |
| `bootstrap_difference_p` | `[ci, p] = bootstrapDifferenceP(La, Lb)` | two bootstrap vectors | 95% CI of the difference, two-sided p | -- |
| `epr_j_models` | `M = eprJModels(J, EPR)` | per-fly J, EPR | linear (b + a|J|) vs quadratic (b + cJ^2) with free intercept, RSS, exponent | G10 |
| `route_duration_observed` | `recs = routeDurationObserved(x, Theta=0.5, LightDark=true)` | 1-s counts | one record per sleep bout: phase, sleep_min, wake_before_min, wake_after_min, q_in, q_out (NaN = not classified) | G11 |
| `route_duration_hidden` | `recs = routeDurationHidden(z, ph, MinSleep=5, MinGap=10)` | one hidden path | same record layout for H3 visits | G11 |
| `route_duration_D` | `D = routeDurationD(recs, relation, phase)` | records | D = mean ln t(quiet) - mean ln t(active), NaN if < 2 of either | G11 |
| `landscape` | `Ls = energyLandscape(P)` | 3 x 3 labelled | E2-E1, E3-E1, branching, barriers, F_SW, F_WS (= xn_pred, yn_pred) | G11 |

`epr_count_bound` is identical to `epr_count_based.m` (rule `'bound'`, the user's existing function). The
two agree to 8.9e-16 over 1e5 random count matrices (ledger). Reuse that file; do not write a second version.

## 3. Conventions (each one has caused or nearly caused an error)

- **Indices.** Python and the golden files are 0-based. MATLAB: add 1 to states (H1 = 1, H2 = 2, H3 = 3), to
  phases, and to every time index or bout start.
  - Symbols: Python Q = 0, W = 1, S = 2, missing = -1, wake of unknown type = -2 (`WAKE_UNKNOWN`). The `.mat`
    export and the suggested MATLAB coding are 1 = Q, 2 = W, 3 = S, with **0 = missing** and **4 = wake of unknown
    type**.
  - A missing minute contributes emission likelihood 1 and is NOT a wake minute in `cObsMinutes`.
  - A wake-of-unknown-type minute contributes P(Q | H) + P(W | H); see "Wake minutes without a Q/W reference"
    below.
- **Time.** Count files start at ZT0, one sample per second.
  - Minute m covers seconds [60m, 60m + 60).
  - Light = ZT 0-12. For minutes the test is `mod(60*m, 86400) < 43200`; for a second s it is
    `mod(s, 86400) >= 43200` -> dark.
  - Sequences are truncated to whole days (1440 min) before fitting.
- **Phase arrays.** `ph(t)` selects the transition matrix for the step INTO minute t; `ph(1)` is never used for a
  transition. Every transition, route and pair count is attributed to the phase (or 4 h window) of its ARRIVAL
  minute.
- **State labels.** After a K = 3 fit:
  - H3 = argmax of P(S | state);
  - H2 = the larger P(W | state) of the other two;
  - H1 = the remaining one.
  - Apply the same permutation to every phase's matrix.
- **Sign.** Positive = forward cycle H3 -> H2 -> H1 -> H3 (Sleep -> Wake -> Reg -> Sleep). This is the SAME cycle
  as A -> Q -> S -> A, and the same convention as the Methods v6 (Eq. 20, 21).
  - `J_cycle = (N32 - N23)/h`.
  - `A = log[P32 P21 P13 / (P31 P12 P23)]` (1-based labels).
  - `c_pred = (Yn - Xn)/2`.
  - For a 3-state stationary chain, sign(c_pred) = sign(A) = sign(J) exactly. The exception is a cycle edge that
    is numerically zero (|J| < 1e-12): there is no cycle and the sign is undefined. Report such blocks as
    "no cycle".
  - Caution: `sqw2_report.py` printed J in the REVERSE direction (H1 -> H2 net). `sqw_ref.model_energetics` and
    everything after the count-J ledger entry use the forward sign.
- **Genotype order (presentation rule).** CS -> per0LD -> per0DD -> fumin -> sss (wild type, clock mutant, sleep
  mutants); light before dark within a genotype. Use it for every axis, legend, table and output loop. Details:
  FIGURE_RECIPES.md, common plotting style.
- **Units.** Steps are minutes. Rates are converted to per hour by x60 (model quantities) or by dividing counts
  by hours in the phase.
- **Non-finite values in golden JSON** are the strings `"Inf"`, `"-Inf"`, `"NaN"`. MATLAB's `jsondecode` returns
  them as char; convert with `str2double`.

## 4. Algorithms: the details a port gets wrong

- **Bridging (lab convention b).**
  - A single active second is reset to inactive ONLY if the inactive runs on BOTH sides are >= 300 s.
  - Flank lengths are measured on the UNBRIDGED trace, so 300+1+300+1+300 becomes one 902 s bout. 301+1+100 is
    not bridged; 300+2+300 is not bridged.
  - A run touching the start of the trace counts as a flank.
  - G1 covers all of these cases.
- **Wake-bout events.**
  - Only wake bouts >= 600 s between two sleep bouts are used.
  - ref = active fraction in the middle 300 s, `[mid-150, mid+150)` with `mid = floor((w0+w1)/2)`.
  - Exit and approach windows are seconds 61-300 (`[w0+60, w0+300)` and `[w1-300, w1-60)`).
  - r_f for a phase = median of ref over events whose waking OR onset is in that phase; it needs >= 3 events.
- **Symbols.**
  - S if >= 30 s of the minute lie in a sleep bout.
  - Otherwise Q if the active fraction < theta * r_f (strict).
  - Otherwise W.
- **Wake minutes without a Q/W reference (RECOMMENDED handling: partial observation).**
  - **When it happens.** A phase whose fly has fewer than 3 qualifying wake bouts (>= 10 min, between two sleep
    bouts) has no reference r_f, so its wake minutes cannot be called Q or W. Seen in 15 fly-phases: fumin 4,
    sss 11, mostly flies awake through most of a phase.
  - **What the analyses so far did.** Coded those minutes as missing (-1, likelihood 1 for every state).
  - **Why that is wrong in principle.** It discards the fact that the minute is NOT sleep. The decoder then places
    H3 there whenever the dynamics allow, which produced spurious H3 visits < 5 min in fumin and sss.
  - **Recommended.** Code them as wake of unknown type (Python -2, MATLAB 4):
    `sqw_symbols(..., unknown_wake=WAKE_UNKNOWN)`. Emission likelihood for state H is P(Q | H) + P(W | H). In the
    EM emission update, each such minute's posterior weight gamma(H) is split between Q and W in proportion to
    P(Q | H) and P(W | H), so the minute counts toward wake but not toward either type.
  - **Implementation.** `emission_lik`, `em_step` and `em_fit` in `sqw_ref.py`. Checked:
    - with no -2 minutes, results are bit-identical to before (golden G3 unchanged);
    - with -2 minutes, the batched `em_fit` equals iterated `em_step` (differences ~1e-15) and the likelihood never
      decreases over 30 EM steps;
    - golden G9 holds the test case.
  - **Downstream.**
    - Runs: for wake runs include -2 in the set, `runs_censored(y, ph, {0, 1, -2}, ...)`; Q-only and W-only runs
      are cut (censored) at a -2 minute.
    - `observable_pair_counts` and `c_obs_minutes` skip -2, as they skip -1.
  - **Not yet applied to the fly results.** Applying it needs refits of the affected flies (12 fumin/sss flies).
    A sensitivity check excluding those fly-phases (PRESPEC 20i) left every result unchanged, so the conclusions
    do not depend on it; a final rerun or the MATLAB version should use it.
- **EM.**
  - Batched Baum-Welch with scaled forward-backward.
  - Transition expected counts are split by phase of arrival; emissions are pooled over phases.
  - Initial values: Dirichlet(1) for pi, every row of P, and every row of B. The first `floor(starts/2)` starts
    use P <- 0.5 P + 0.5 I.
  - Stop when, for EVERY start, (ll - ll_old) < 1e-7 * |ll|, or after 1000 iterations; keep the best final
    log-likelihood.
  - A port cannot reproduce numpy's random starts. Validate `hmmEmStep` exactly against G3. On real data,
    check that the port reaches the same optimum: same log-likelihood to ~1e-4 for most flies. A port finding
    a HIGHER likelihood is legitimate, not an error.
- **Do not use `hmmtrain` / `hmmviterbi`** (Statistics Toolbox). They cannot share emissions across phases or
  switch transition matrices by phase. Icode's bout suite deliberately avoids toolboxes; keep that.
- **Model order.**
  - k = (K-1) + nP K (K-1) + 2K.
  - BIC with n = number of minutes (primary) and n_eff = number of symbol changes (sensitivity).
  - Flag a fly if K = 2 or K = 4 beats K = 3 by > 10.
  - Routing is defined only at K = 3.
- **Degeneracy.** H1 is not identified if its emission row is within total variation 0.05 of H2's or H3's.
  Exclude such flies from routing; report how many.
- **FFBS.**
  - Forward pass first.
  - Sample the last state from the final alpha.
  - Then backward: p(z_t = i | z_{t+1} = j) is proportional to alpha_t(i) P_{ph(t+1)}(i, j).
  - 200 paths per fly.
- **Viterbi** is kept only as a comparison. Ties go to the lowest index.
- **Routes** are matched on the run-compressed path:
  - H3 H1 H2 -> N312;
  - H3 H2 -> N32;
  - H2 H1 H3 -> N213;
  - H2 H3 -> N23.
  - Excursions H3 H1 H3 and H2 H1 H2 are not routes.
  - Pair counts N_ij are single-step changes on the same compressed path.
- **Jensen.** Compute nonlinear statistics (Xn, Yn, c, EPR) on EACH posterior path, then average. Do not compute
  them on the path-averaged counts. J_cycle is linear and is unaffected.
  - Per path, Xn needs >= 3 routes in its denominator.
  - A fly's per-path mean is reported only if the statistic is defined on >= 50% of its paths.
- **Rule of three.** A pair seen in one direction only, count k:
  - the term is (k-3) ln(k/3) if k > 3, else 0;
  - equivalently, b* = min(3, k) enters BOTH the prefactor and the log;
  - never substitute a flat 3;
  - assert EPR >= 0.
- **Stationary distribution.** `sqw_ref.stationary` uses the eigenvector for eigenvalue 1. The Methods v6
  QR least-squares solve of [P' - I; 1'] pi = [0; 1] is better conditioned for near-reducible P; prefer it in
  MATLAB. Checked: the two agree to <= 2.4e-15 on G4 and <= 6.6e-14 on every fitted 20b matrix.
- **Model energetics.** An edge with P < 1e-10 in exactly one direction makes EPR and A infinite. Report the
  infinite value; do not floor it.

## 5. Testing the port

- One `test_*.m` per group of functions, run by name from the folder (Icode CLAUDE.md, Testing). Read the golden
  JSON with `jsondecode(fileread(...))`.
- **Exact (deterministic) tests:**
  - `sleepBouts` (G1);
  - `sqwSymbols` (G2: y, refs, bouts);
  - `hmmForward` log-likelihood and `hmmEmStep` (G3, tolerance 1e-10 relative);
  - `hmmViterbiPhased` (G3);
  - labelling, degeneracy, routing, energetics (G4, 1e-12);
  - counting and the EPR bound (G5, exact or 1e-12).
  - observable pair counts and round trips (G7, exact).
  - dwell-time runs, Kaplan-Meier, hazard (G8, exact).
  - Onsager fit and EPR-J models (G10, 1e-12); bootstrap statistically.
  - route-duration records, D and the landscape (G11, exact).
- **Statistical tests:**
  - `hmmFFBS`: posterior mean route counts over 4000 paths within 4 SE of G6;
  - `hmmFitPhased`: on G3's inputs, the final log-likelihood is >= that of the G3 parameters after many
    `hmmEmStep` iterations;
  - `hmmSimulatePhased`: empirical transition and emission frequencies within sampling error.
- **End-to-end:** on one real fly, reproduce the ledger's 20b values from the stored fitted parameters
  (`sqw_out/fit2_0.5.jsonl`, not in the repo; export it with `sqw2.py mat` or ask for it).
- **Do not "validate" with seeds.** The same seed gives different draws in MATLAB and Python (Icode CLAUDE.md).

## 6. What changed relative to the previous MATLAB pipeline (Methods v6, `compute_HMM_energetics.m`)

| | previous | now | why |
|---|---|---|---|
| Observable | binary S/W per second | S/Q/W per minute | H1 had no distinct emission, so routing was not identified (PRESPEC 1-2) |
| Fit unit | each 4 h window | each fly, whole record, light and dark matrices, shared emissions | one window cannot resolve routing (PRESPEC 20, V1 failed) |
| Model order | AICc, m = 2..7, per window, windows with m != 3 dropped | BIC, K = 2..4, per fly, flagged not dropped | per-fly power; 3 states chosen in 76-96% of flies |
| Fitter | `hmmtrain` | own phase-split Baum-Welch | shared emissions, phase-specific P |
| Decoding | Viterbi per window, windows stitched | FFBS over the whole record, 200 paths; Viterbi for comparison | Viterbi biased when emissions are degenerate; no seams |
| Zero counts | P_min 'drop' (+Inf -> 0) | rule-of-three lower bound | `epr_count_based.m` header |
| Window extension | windows extended to finish a bout | never | length bias (PRESPEC 8b) |

The per0 two-subgroup result (Methods v6, "Detection of routing subgroups within per0") was derived from
routing in the binary model. There, a fly's side of the diagonal is predicted by how often its degenerate H1
copied the sleep emission row (Spearman -0.94; ledger, post hoc). It does not reappear in the S/Q/W model.

## 7. Files

- Code to port: `sqw_ref.py`.
- Its check: `verify_ref.py` (needs the scratchpad data; it is a record, not something to rerun).
- Tests to satisfy: `golden/`, made by `make_golden.py`.
- Raw counts and `.mat` files are NOT in the repository; use local copies. Layout expected by the scripts:
  `counts/<group>/**/*.txt`, with groups CS, fumin, per0LD, per0DD, sss.
- Analysis scripts, for provenance only (see `SCRIPTS_INDEX.md`): `sqw2.py` (20b), `sqw3.py` (20c),
  `sqw4.py` (20d), `sqw5.py` (20e), `sqw2_report.py`, `sqw2_countJ.py`.

## 8. Per-fly round-trip table

`roundtrip_table.py [theta]` writes `roundtrip_per_fly_theta<theta>.csv` using only `sqw_ref.py` functions and the
stored 20b fits: one row per fly x phase x level (observed, hidden_all, hidden_matched). Columns: genotype,
fly_index, file, phase, level, theta, n_round_trips, frac_* (4 types), rate_per_h_* (4), p_quiet_entry,
p_quiet_exit, dependence_* (4), frac_paths_used. It reproduces the 20g / 20g-B per-fly values to 5e-5 (rounding of
the CSV). A MATLAB port should reproduce this table from the same fitted parameters.

## 9. Per-fly D: does the route depend on bout duration? (PRESPEC 20k; Fig. 13d)

Inputs per fly: the 1-s counts `x` (observed level), and its fitted 20b model with 200 posterior paths (hidden
level). theta = 0.5.

### 9.1 Records
**Observed.** `routeDurationObserved`. Take every bridged sleep bout k that has a sleep bout before and after it.
1. Durations in minutes:
   - sleep_min = length of bout k;
   - wake_before_min = start of k minus end of k-1;
   - wake_after_min = start of k+1 minus end of k.
2. Phase = phase of the onset (light if mod(onset_s, 86400) < 43200; per0DD a single phase).
3. q_in:
   - only if wake_before >= 600 s;
   - TRUE if the active fraction in [onset - 300 s, onset - 60 s) < theta * r_f(phase).
4. q_out:
   - only if wake_after >= 600 s;
   - TRUE if the active fraction in [end + 60 s, end + 300 s) < theta * r_f(phase).
5. Otherwise the field is NaN (not classified).
6. Skip the fly-phase if r_f is undefined.

**Hidden.** `routeDurationHidden`, on each path z.
1. Use every H3 run >= 5 min that has a previous and a next H3 run.
2. gap_before = start of this run - end of the previous H3 run; gap_after = start of the next H3 run - end of this
   one.
3. q_in = (state before == H1) if gap_before >= 10, else NaN; q_out = (state after == H1) if gap_after >= 10,
   else NaN.
4. Phase = phase of the run's first minute.

### 9.2 Relations
| relation | route | duration |
|---|---|---|
| R1 | q_out | sleep_min |
| R2 | q_in | wake_before_min |
| R3 | q_in | sleep_min |

### 9.3 D for one fly and phase (`routeDurationD`)
1. Keep the records of that phase whose route is classified.
2. Split them into quiet (route TRUE) and active (FALSE).
3. D = mean(log(t_quiet)) - mean(log(t_active)), natural log, t in minutes.
4. D is NaN if either group has fewer than 2 bouts.
5. Hidden level: compute D on each of the 200 paths. The fly's D = mean over paths where it is defined, and is
   NaN if it is defined on fewer than 50% of paths.
6. Interpretation: D < 0 means quiet routes go with shorter bouts. exp(D) is the ratio of geometric-mean durations,
   e.g. D = -0.29 -> 0.75.

### 9.4 Group statistic and nulls
- Group value = plain mean of D over flies in the genotype x phase (each fly once).
- **Observed null.**
  - For each fly separately, randomly permute its route labels among its own classified bouts in that phase,
    recompute its D, and average over flies. Repeat 2000 times.
  - The permutation keeps each fly's durations and its numbers of quiet and active routes; only the pairing is
    broken.
  - p = 2 * min(fraction of null >= real, fraction <= real).
- **Hidden null.**
  - For each fly: simulate a record from its fitted model (same length and phases), decode 20 FFBS paths with the
    true parameters, and compute its path-averaged D.
  - Repeat 400 times. The r-th group null value = mean over flies of their r-th simulated D.
  - p = 2 * min((1 + #null >= real)/(1 + 400), (1 + #null <= real)/(1 + 400)).
  - Do not use only 40 replicates: the smallest attainable p is then 0.049, and nothing survives Holm (PRESPEC 20k
    amendment A).
- Holm correction over the 9 genotype x phase blocks, separately per relation and per level.

### 9.5 Figure (`fig_perfly_D.py`; FIGURE_RECIPES I6b)
- One point per fly (D), jittered; black bar = group mean.
- Grey box = 2.5-97.5 percentiles of the group null.
- * for Holm p < 0.05.
- Rows observed / hidden; columns R1 / R2 / R3.

### 9.6 Checks
- CS dark R1: observed mean D -0.293 (79% of flies < 0, Holm p < 0.001); hidden -0.275 (76%, Holm p 0.045).
- per0LD light R1 observed: -0.228 (Holm 0.04).
- Exact unit test: `golden/G11_route_duration.json`, which covers hidden records on a fixed path, observed records
  from the G2 trace, and D for a hand-checkable set of 7 bouts (R1 = -1.3088).
