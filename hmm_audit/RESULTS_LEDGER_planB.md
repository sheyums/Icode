# RESULTS LEDGER — Plan B feasibility (branch planB-feasibility; append-only)
Separate from the shared RESULTS_LEDGER.md. Every entry gives the PRESPEC commit, the code commit and the date.

## PRESPEC 20t F1 + F5 — synthetic recovery of clock / homeostat modulation and of the EP split (2026-11-01; PRESPEC 3465f69, code 1041661, report script bfcae17; output planB_f1f5.json, planB_f1f5_report.log)
Narrowest block: one synthetic fly's whole record (4 d DD-type; 2 d LD2-CH; 3 d LD3-step-C); 25 flies per scenario.
Sample-size-matched null: N0 (4 d) and N0-2d (2 d), same fitting procedure. Run time 4153 s on 4 cores. The first launch
(2026-10-31) was killed when the session paused, before writing any output; this is a re-run of the same commit.
```
Detection rates (LRT > nominal chi2 5%): clock rows Q/W/S | homeostat rows Q/W/S
N0          n 25 | clock 0.08 0.04 0.12 | homeo 0.16 0.24 0.08
DD-CH       n 25 | clock 0.24 0.04 0.32 | homeo 0.24 0.24 0.44
DD-CH-half  n 25 | clock 0.12 0.08 0.08 | homeo 0.24 0.24 0.32
DD-C        n 25 | clock 0.28 0.00 0.48 | homeo 0.04 0.04 0.16
DD-H        n 25 | clock 0.00 0.08 0.12 | homeo 0.20 0.16 0.36
LD2-CH      n 25 | clock 0.16 0.04 0.28 | homeo 0.08 0.24 0.24
LD3-step-C  n 25 | clock 0.44 0.00 0.52 | homeo 0.04 0.08 0.12
N0-2d       n 25 | clock 0.08 0.00 0.20 | homeo 0.12 0.08 0.20
  N0 95th pct LR: clock 10.39 9.07 10.64 | homeo 8.20 7.73 7.30   (nominal 9.49 / 5.99)
  N0-2d 95th pct LR: clock 12.27 6.89 12.57 | homeo 9.56 6.66 15.41
Estimates (median over flies): 02 = Q->S, 21 = S->W
DD-CH  02 Eb ratio 1.16 [IQR 0.10,1.99] sign ok 0.76 | a ratio 0.99, phase err 1.7 h
DD-CH  21 Eb ratio 1.69 [0.83,2.64] sign ok 0.96     | a ratio 1.06, phase err 1.3 h
LD2-CH 02 Eb ratio 1.98 [1.14,3.95] sign ok 0.84     | a ratio 1.74, phase err 2.0 h
LD2-CH 21 Eb ratio 1.82 [0.85,3.24] sign ok 0.92     | a ratio 1.72, phase err 1.9 h
N0 (truth a = 0): a_hat median 0.36 (Q->S), 0.24 (S->W); N0-2d 0.44, 0.35
tau_r, tau_d (truth ~6 h, ~3 h): profiled medians 16-32 h and 4-32 h in every scenario (not recovered)
r(H, clock) multiple correlation, true H: DD-CH 0.47, DD-C 0.43, DD-H 0.22, LD2-CH 0.46, N0 0.30
F5 (nats/h, medians; sigma_na SIGN CORRECTED, see erratum)
N0      hk true 0.091 hat 0.146 ratio 2.19 | na true 0       hat +0.0018
DD-CH   hk true 0.128 hat 0.211 ratio 1.64 | na true +0.0012 hat +0.0065 ratio 6.08 [3.42,8.40] | na share of total 1.2%
DD-CH-h hk true 0.142 hat 0.306 ratio 1.69 | na true +0.0003 hat +0.0028 ratio 7.92
LD2-CH  hk true 0.170 hat 0.393 ratio 2.46 | na true +0.0010 hat +0.0131 ratio 12.97 [4.19,27.50] | na share 0.5%
noise floor: N0 fitted na / DD-CH true na = 1.54
```
VERDICTS by the pre-specified criteria:
- **F1 FAIL.** (1) detection of truly modulated terms: min 0.24 in DD-CH (best 0.44), 0.08 in LD2-CH; required >= 0.80.
  (2) false positives: homeostat 0.24 in the unmodulated W row and up to 0.24 in N0 (the profiled time constants make the
  chi-square anti-conservative, as anticipated); required <= 0.10. (3) cross-talk: homeostat 0.16 in the S row of DD-C;
  required <= 0.10. (4) E_b ratio 1.16 / 1.69 (DD-CH) and 1.98 / 1.82 (LD2-CH), required [0.75, 1.25]; sign 0.76-0.96.
- **F5 FAIL.** Housekeeping over-estimated x1.64 (DD-CH), x2.46 (LD2-CH); driving over-estimated x6-13; the driving
  term's noise floor (N0) is 1.5x the true DD-CH value. The true driving share is only 0.5-1.2% of the total.
- F6 NOT RUN, as pre-specified (it runs only after F1 passes).

Reading (no new test):
- The cause is the number of events, not the fitting. Per fly and row, 4 days give roughly 50-250 transitions out of Q
  or S. Under the null the fitted amplitude is already 0.24-0.44, as big as the effect sought (0.4). Per-fly
  identification of WHICH transition a regulator acts on is out of reach at these record lengths.
- Group-level recovery is better: in DD-CH the median clock amplitude ratio is 0.99 / 1.06 and the phase error
  1.3-1.7 h. This is an observation from the same run, NOT a pass criterion. Reading per-fly estimates as a group
  (pooled or hierarchical) would be a NEW design and needs its own PRESPEC.
- The process-S time constants are not identified (profiled at 16-32 h against a truth of 6 h / 3 h). b is therefore
  confounded with its time scale; E_b = b * SD(H) was chosen to be scale-free, and it is still biased high.
- Plug-in EP is biased upward (a positive quantity estimated from noisy rates; even N0 is x2.2). The driving part is
  about 1% of the total, so a measured "cost of being driven" would be below the noise floor at these data sizes.
- LD3-step-C (diagnostic): light masking fitted with cos + H gives spurious homeostat detection 0.04 / 0.12 (Q / S),
  similar to the null level. The clock term absorbs most of the step (a_hat ratio 1.6-1.8).

ERRATUM (found while reading the output, before any decision): the code's sigma_na had its log argument inverted,
and so did the formula in PRESPEC_20t.md line 62. The correct Esposito & Van den Broeck term is
sum p_i k_ij ln[p_i pi_j / (p_j pi_i)], the one the PRESPEC's own "total = hk + na" requires. The sign flips truth and
fit identically, so every ratio and verdict above is unchanged; only the printed sign of na was wrong (shown corrected
above). The code line was fixed in the commit that adds this entry; no re-run was needed.

## PRESPEC 20t F2 — separability of clock and homeostat, real flies (2026-11-01; PRESPEC 3465f69, code 6698097; planB_f2.json, planB_f2.log)
Narrowest block: one fly's whole record (LD 3-4 d, per0 DD 3-4 d). Data: the user's upload of 2026-11-01 (CS 25, per0LD 27,
fumin 53, sss 57, per0DD 25 files). Identity check: CS c_obs medians light +0.042, dark +0.147 = 20s-A exactly.
CS DD was NOT in the upload, so F2 for CS DD, the genotype the criterion matters most for, is not done.
```
CS      n 25 | r(6h,3h) median 0.47 IQR [0.33,0.65] max 0.89 | r<0.8 in 0.96 -> PASS | worst grid point: r<0.8 in 0.88
per0LD  n 27 | r(6h,3h) median 0.59 IQR [0.41,0.70] max 0.87 | r<0.8 in 0.89 -> PASS | worst grid point: r<0.8 in 0.85
fumin   n 53 | r(6h,3h) median 0.58 IQR [0.41,0.77] max 0.96 | r<0.8 in 0.81 -> PASS | worst grid point: r<0.8 in 0.74
sss     n 57 | r(6h,3h) median 0.59 IQR [0.42,0.73] max 0.89 | r<0.8 in 0.91 -> PASS | worst grid point: r<0.8 in 0.82
per0DD  n 25 | r(6h,3h) median 0.28 IQR [0.22,0.34] max 0.52 | r<0.8 in 1.00 -> PASS | worst grid point: r<0.8 in 1.00
```
- F2 PASS in all 5 available genotypes. Real-fly collinearity (0.47-0.59 in LD) is close to the synthetic DD-CH flies'
  (0.47), so F1's low per-fly power was NOT caused by an unrealistically collinear synthetic design.
- Does not change the step-4 decision (F1 failed). In LD the cos term also stands in for the light step, so r here
  measures H vs 24-h periodicity, not H vs the clock alone.
- F3 (needs CS DD) and F4 (needs deprivation recordings): still not runnable; neither was in this upload.

## PRESPEC 20t F3 + F4 — clock-precision bound; deprivation data (2026-11-01; PRESPEC 3465f69, code 7cc22e5 + planB_f4_lost.py in this commit; planB_f3f4.log, planB_f4_lost.log, planB_f3f4.json)
Data: CS_DD_activfiles.zip (25 files, 4.000 d each, the cut DD files) and CS_9hrdeprExp.zip (26 files, 9.513 d each),
both uploaded 2026-11-01; per0 DD as before.
```
F3 CS DD   [rhythmic] flies 25 intervals 75 | tau 24.03 h | s 1.60 h [1.26, 1.92] | B = 2 tau/s^2 = 18.7 k_B/h [13.0, 30.3] | CI ratio 2.32
F3 CS DD   rhythmic 25 of 25 | periods median 23.90 h
F3 per0 DD [all     ] flies 25 intervals 66 | tau 23.57 h | s 7.16 h [6.30, 7.91] | B = 0.9 k_B/h [0.7, 1.2]
F3 per0 DD rhythmic 23 of 25 (!)
F3 criterion: rhythmic 25 >= 15 and CI ratio 2.32 < 4 -> PASS
F4 n 26 | baseline ZT15-24 sleep median 367 min | deprivation night 51 min (range 0-316) | sleep lost median 315 min,
   P10 172, P90 375, range 84-502 | P90/P10 2.18 | CV 0.33 -> PASS
```
- **F3 PASS by the rule, but its rhythmicity half does not discriminate.** The uncorrected chi-square periodogram
  (81 trial periods, autocorrelated 1-min counts) calls 23/25 arrhythmic per0 DD flies "rhythmic". That half of the
  criterion is uninformative and must not be quoted as a rhythmicity result.
- The informative number is the jitter. CS DD s = 1.60 h [1.26, 1.92] against per0 DD 7.16 h. 7.16 h is close to
  the ~6.9 h SD of a uniformly random marker on the wrapped +/-12 h difference, so per0's "bound" of 0.9 k_B/h is the
  method's floor, not a measurement.
- CS DD minimum clock dissipation B = 2 tau / s^2 = 18.7 k_B/h (95% CI 13-30), about 20x the arrhythmic floor and
  about 200-350x the observed S/Q/W cycle's EPR (0.053 k_B/h, 20s-A).
- The centroid marker adds noise, so s is an OVERestimate and B a valid but weak lower bound. Only 3 intervals per fly
  (4-day records); longer records would tighten it.
- **F4 PASS.** The deprivation sits on DAY INDEX 3, ZT15-24 (group sleep 2-12 min/h against 38-50 on days 0-2), i.e.
  THREE baseline days, not two as remembered. 6 recovery days follow (day 9 partial).
  - Located from the group hourly profile before any per-fly number was computed. The window agrees with the user's
    "last 9 h of the night".
  - Rebound is visible but modest: day 4 ZT1-6 43-55 min/h vs 31-51 on baseline days.
- **Caveat for any dose-response:** the spread in sleep lost comes partly from deprivation failing in some flies (one
  slept 316 of 540 min). Dose is then partly a fly property (resistance to arousal), which can confound "dose" with
  sleep drive. A dose-response PRESPEC must address it, e.g. with deprivation efficacy as a covariate, or by
  restricting to flies with >= 80% of the window awake.

## PRESPEC 20u Stage 0 — gates (2026-11-01; PRESPEC ccae0e8 + amendment A 3f3f2cd; code ae2c364; planB_20u_stage0.log/.json)
Go-ahead for tests 1-3 (20u, 20v, 20w) given by Sheyum, relayed by cClaude, 2026-11-01.
```
G1: port wake 4017 (ref 4017), sleep 4022 (ref 4022); flies identical 25/25 -> PASS
     longest inactive run on real traces 11881 s = 3.30 h (amendment A: must be < 12 h)
G2: RiseRatio 12 bins 1.4155 16 bins 1.4022 20 bins 1.3855 24 bins 1.4347 (target 1.416, 1.402, 1.385, 1.435) -> PASS
```
- The Python extractor (sqw_ref.sleep_bouts + edge sleep bouts dropped) reproduces lClaude's MATLAB per0 DD bouts EXACTLY
  per fly (sorted wake and sleep bouts identical, 25/25 flies). The RiseRatio port reproduces 20r-B's observed values.
  Fallback (a) is not needed. Proceeding to stage 1.

## PRESPEC 20w — clock face, descriptive (2026-11-01; PRESPEC ccae0e8 + 20u amendment A item 4; code e6fc5f3; planB_20w.log/.json)
Narrowest block: one 2-h clock bin pooled over the flies and days of a genotype. 2000 fly-bootstrap replicates.
Pooled-chain J (cycles/h) by 2-h bin, point estimates (full table with bands, A, EPR, E3-E1 and c_obs in planB_20w.log):
```
bin h     0-2   2-4   4-6   6-8  8-10 10-12 12-14 14-16 16-18 18-20 20-22 22-24 | range [95% band]
CS LD   +.025 +.017 +.001 +.079 +.115 +.008 +.228 +.146 +.250 +.086 +.088 +.177 | 0.250 [0.220, 0.448]
per0 LD +.150 +.113 +.084 +.096 +.033 +.081 +.040 +.009 +.166 +.122 +.191 +.209 | 0.200 [0.156, 0.369]
per0 DD +.183 +.108 +.190 +.141 +.152 +.194 +.102 +.167 +.193 +.173 +.137 +.162 | 0.092 [0.120, 0.355]
CS DD   +.251 +.077 +.084 +.099 +.019 +.052 +.074 +.063 +.089 +.186 +.086 +.166 | 0.232 [0.172, 0.385]
EPR range: CS 0.251 [0.163, 0.532]; per0LD 0.131 [0.082, 0.333]; per0DD 0.085 [0.094, 0.391]; CSDD 0.231 [0.128, 0.440]
```
- PRE-SPECIFIED READING: no genotype's range band lies entirely above per0 DD's, for J or for EPR, so none is called
  "rhythmic" by this summary.
- The summary is insensitive, and this is visible in the output: the bootstrapped max-minus-min range is inflated by
  noise. per0 DD's band [0.120, 0.355] lies entirely ABOVE its own point estimate 0.092, which a 12-bin range at this
  bin width does to any genotype. Stated as a limitation of the chosen summary, not fixed (no re-analysis).
- Descriptively the profiles agree with C3: CS LD J is ~0 at ZT0-6 and 10-12, peaks at ZT12-14 and 16-18
  (+0.23, +0.25), and is forward through the night. CS DD is highest at CT0-2 and 18-24 (late subjective night).
  per0 DD is flat within its bands (0.10-0.19 in every bin). per0 LD is highest late night and early morning.
- E3-E1 stays negative in almost every bin of every genotype (sleep is the lowest-energy state round the clock);
  CS DD's minimum magnitude is at CT8-12. Single-bin values have wide bands; no bin-wise test was planned or run.
- CS DD bins are clock time from the cut start; free-running drift (<= ~2 h over 4 days) smears them (stated in 20w).

## PRESPEC 20v — waiting-time irreversibility beyond the observed chain (2026-11-01; PRESPEC ccae0e8 + 20u amendment A item 4; code e6fc5f3; planB_20v.log/.json)
Narrowest block: one fly's whole DD record (LD: one fly's light or dark minutes pooled over days). Null A = 20p's
(own mc_fit chain, 200 records, missing reinserted). Excess E_WT = sigma_WT(obs) - median(null), nats/h.
```
per0DD  whole n 25 [PRIMARY] | sigma_WT obs median 0.2756 | E_WT median +0.0827 | Wilcoxon p 0.000556 Holm 0.000556 | flies p<.05 0.28 (binom 0.000169) -> SUPPORTED
CSDD    whole n 25 [PRIMARY] | sigma_WT obs median 0.2604 | E_WT median +0.0782 | Wilcoxon p 8.8e-05 Holm 0.000176 | flies p<.05 0.32 (binom 1.96e-05) -> SUPPORTED
CS      light n 25 [descr.] | 12-h null: E +0.0281 p 0.442 frac 0.08 | 6-h null: E +0.0345 p 0.578
CS      dark  n 25 [descr.] | 12-h null: E +0.1367 p 0.0013 frac 0.04 | 6-h null: E +0.1284 p 0.00309
per0LD  light n 27 [descr.] | 12-h null: E +0.0269 p 0.0463 frac 0.15 | 6-h null: E +0.0293 p 0.0655
per0LD  dark  n 27 [descr.] | 12-h null: E +0.0295 p 0.269 frac 0.04  | 6-h null: E +0.0379 p 0.336
fumin   light n 49 [descr.] | 12-h null: E +0.0323 p 0.01 frac 0.02   | 6-h null: E +0.0399 p 0.0131
fumin   dark  n 53 [descr.] | 12-h null: E +0.0530 p 6.4e-05 frac 0.11 | 6-h null: E +0.0611 p 8.0e-05
sss     light n 50 [descr.] | 12-h null: E +0.0036 p 0.258 frac 0.02  | 6-h null: E +0.0115 p 0.192
sss     dark  n 53 [descr.] | 12-h null: E -0.0000 p 0.449 frac 0.02  | 6-h null: E -0.0045 p 0.374
```
- VERDICT by the pre-specified rule: SUPPORTED in both primary blocks (per0 DD, CS DD). Dwell times in Q and W depend
  on the direction of passage (W->Q->S vs S->Q->W; Q->W->S vs S->W->Q) beyond what each fly's own first-order chain
  produces. The excess is ~0.08 nats/h, comparable to the observed chain's own EPR (per0 DD sigma_MC median 0.111;
  CS DD 0.053; 20s-A).
- Agreement with 20p: per0 DD carries excess irreversibility again, by an independent statistic (dwell shape
  rather than block frequencies). CS DD also qualifies here. In 20p, CS DD's whole-record row had qualified but its
  designated 12-h halves had not.
- MAJOR CAVEAT, NOT EXCLUDED by this design: the mixed pairs are conditioned on transitions into and out of S, and S
  is scored by the 5-min criterion. The wake minute(s) before a sleep onset contain the start of the inactivity that
  becomes sleep, while waking is followed by activity. That scoring asymmetry can make Q/W dwells next to S
  direction-dependent, and a first-order chain on the symbols does not reproduce it. 20p excluded this for per0 DD
  with wake-only blocks; 20v has no wake-only guard (cClaude flagged the same issue for S dwells; it applies to the
  conditioning too). Until a guard is run, the 20v excess must be worded "consistent with hidden irreversibility; a
  sleep-scoring origin is not excluded". The natural guard (not run, not pre-registered) is the same statistic on
  records simulated at 1 s and scored through sqw_symbols, i.e. a scoring-matched null.
- sss (descriptive): E_WT ~0 in light and dark, so no waiting-time evidence for a hidden cycle without a visible current.
- LD (descriptive): excess in CS dark and fumin light and dark, not CS light, per0LD or sss. It is unchanged by the
  6-h-rate null, so within-phase drift does not explain it.
- Not a bound by itself: the lower-bound statement is sigma_MC + E_WT per fly, which was not computed here.

## PRESPEC 20w-B — clock-face cosine test (2026-11-01; PRESPEC 43d0e34, code 5bdf525; planB_20wB.log/.json)
Narrowest block: each genotype's whole recording, pooled in 2-h clock-time bins (20w's data, bins and pooling).
Null: 999 per-fly circular rotations by uniform whole-minute offsets in [0, 24 h); sample-size matched by construction.
```
check: max |J_vec - model_energetics| = 8.33e-15
CALIBRATION: fraction p < 0.05 = 0.045 over 200 per0-DD-matched constant-chain datasets (gate [0.02, 0.09]) -> PASS
POWER at per0 DD size vs CS LD's pooled modulation (no heterogeneity, optimistic): 0.80 (100 datasets)
CS      n 25 | J amplitude 0.0774 (null 95th 0.0618) peak 16.2 h | p 0.005 | Holm 0.015 -> CLOCK-LOCKED MODULATION | c_obs amp 0.0691 (null 95th 0.0517) peak 17.5 h p 0.003
CSDD    n 25 | J amplitude 0.0628 (null 95th 0.0551) peak 22.9 h | p 0.015 | Holm 0.030 -> CLOCK-LOCKED MODULATION | c_obs amp 0.0468 (null 95th 0.0438) peak 23.4 h p 0.034
per0LD  n 27 | J amplitude 0.0711 (null 95th 0.0725) peak 22.7 h | p 0.059 | Holm 0.059 -> not detected | c_obs amp 0.0422 (null 95th 0.0516) peak 23.8 h p 0.140
per0DD  n 25 | J amplitude 0.0058 (null 95th 0.0503) peak 18.3 h | p 0.959 | NEGATIVE CONTROL | c_obs amp 0.0042 (null 95th 0.0312) p 0.953
fumin   n 53 | J amplitude 0.0700 (null 95th 0.0455) peak 11.0 h | p 0.001 | descriptive | c_obs amp 0.0515 (null 95th 0.0564) peak 13.7 h p 0.091
sss     n 57 | J amplitude 0.0334 (null 95th 0.0242) peak 16.5 h | p 0.002 | descriptive | c_obs amp 0.0893 (null 95th 0.0604) peak 17.5 h p 0.002
per0 DD null (p >= 0.05): no flag
```
- Calibration gate PASSED (0.045, close to nominal). The per0 DD control is clean (p 0.96; amplitude 0.006 against a
  null 95th of 0.050), so there is NO flag.
- VERDICT (Holm over CS LD, CS DD, per0 LD): clock-locked modulation of the cycle flux J in CS LD (Holm 0.015, peak
  ZT16) and CS DD (Holm 0.030, peak CT23); per0 LD not detected (Holm 0.059; amplitude 0.071 against a null 95th of
  0.0725, a near miss, not a positive).
- The CS LD positive is light plus clock; this test does not separate them. CS DD has no light, so its positive is the
  clock's own modulation of J (flies sharing phase after release into DD), DILUTED by the free-running smear. Its
  amplitude (0.063) is therefore a lower estimate.
- per0 LD's near miss is consistent with light-driven modulation without a clock, but it is not a result.
  Optimistic power is 0.80 against CS LD's amplitude at per0 DD's size; per0 LD's own amplitude (0.071) is similar.
- Secondary c_obs (outside the family) agrees in sign and phase: CS 0.003, CS DD 0.034, per0 LD 0.14, per0 DD 0.95.
- Descriptive: fumin (p 0.001, peak ~ZT11) and sss (p 0.002, peak ~ZT16-17) also have a clock-locked J modulation,
  sss around a mean flux of ~0.
- This supersedes the 20w range summary as the clock-face statistic. 20w stands as recorded (descriptive, no
  genotype rhythmic by its rule).

## Note (2026-11-01): 20u stage 1 crashed; numerical fix, no result seen
The first 20u run (code 0fd7221) died in stage 1 with ZeroDivisionError inside the forward-backward kernel
(planB_20u_crash1.log). Cause: a random EM start can give an inactive block whose largest eigenvalue, raised to a run
length up to 11,881 s, underflows to 0, so a normalisation was 0. No fit, adequacy number or null had been produced
or printed. Fix (purely numerical, same model, starts, orders and rules):
- lambda_max^(L-1) is factored out of each run's propagation, with (L-1) log lambda_max added to the log-likelihood;
- the counts use the identically rescaled G and Z;
- a zero normalisation now returns logL = -Inf, ending that start.
Re-validated against the brute-force per-second forward-backward: logL agrees to 3e-14, expected counts to 1.5e-11,
and recovery is unchanged (planB_20u_validate.py). A fast-only inactive block facing a 12,000-s silence now gives a
finite logL with counts summing to T - 1. Re-run with the fixed code (commit below).

## PRESPEC 20v amendment F — F1 gate: time symmetry of the scoring (2026-11-01; amendment 87d9d77, code 8d5192a; planB_20vF_f1.log/.json)
```
F1 per0DD  flies 25 | minutes differing total 0 of 131040 | worst fly agreement 1.00000 -> PASS
F1 CSDD    flies 25 | minutes differing total 0 of 144000 | worst fly agreement 1.00000 -> PASS
```
- reverse(sqw_symbols(x)) == sqw_symbols(reverse(x)) in EVERY minute of every per0 DD and CS DD whole-day retained trace
  (light_dark=False, theta 0.5). The one-second r_f rounding at mid = (w0+w1)//2 changed no symbol here. The scoring
  is exactly time-symmetric on these data, so F0's premise holds in the code. F2 follows.

## PRESPEC 20x — sleep deprivation, within-fly contrasts at ZT0-6 (2026-11-01; PRESPEC 1502a20, code 60fce8a; planB_20x.log/.json; post hoc lines planB_20x_posthoc.log)
Narrowest block: one fly's ZT0-6 on day 4 (360 min). Matched null: placebo day 2 vs days 0-1 (same flies, same window).
r_f from baseline days 0-2 only. 26 flies, all 9.5127 d; one fly without a dark-phase baseline reference.
```
STAGE 0: sleep in ZT0-6, day 4 - days 0-2: median +42.0 min [95% CI +21.2, +57.1] (baseline median 229 min) | Wilcoxon p 6.0e-05 -> PASS
c_obs              n 26 | baseline -0.027 | Delta -0.062 [-0.151, +0.004] p 0.084 Holm 0.168 | Delta_pl -0.014, paired p 0.635 | trend p 0.745 -> NOT SUPPORTED
count J (cyc/h)    n 26 | baseline -0.056 | Delta -0.056 [-0.194, +0.056] p 0.246 Holm 0.246 | Delta_pl +0.000, paired p 0.957 | trend p 0.626 -> NOT SUPPORTED
count EPR (nats/h) n 26 | baseline +0.044 | Delta +0.037 [-0.003, +0.152] p 0.024 Holm 0.071 | Delta_pl +0.091, paired p 0.617 | trend p 0.059 -> NOT SUPPORTED
no-cycle / rule-of-three kinds (condition:kind -> flies): R:i 3, R:ii 2, P2:i 1, P2:iii 1, P01:ii 1
transitions, pooled chain, recovery/baseline [fly-bootstrap 95% CI], Holm over 6:
  Q->W 0.94 [0.83, 1.05] 0.542 | Q->S 1.16 [0.98, 1.35] 0.243 | W->Q 1.17 [1.05, 1.30] 0.020 | W->S 1.49 [1.24, 1.80] <0.001
  S->Q 0.90 [0.71, 1.12] 0.542 | S->W 0.75 [0.58, 0.94] 0.040
dose association (Spearman, sleep lost vs Delta): c_obs +0.00 (all) / +0.31 (efficacious, n 19); J +0.03 / +0.32; EPR +0.19 / -0.31; all p >= 0.19
descriptive median Delta vs same window of days 0-2 (sleep min | c_obs | J | EPR):
  d4 ZT6-12 +3.7 | +0.010 | +0.028 | +0.188     d4 ZT12-24 +33.9 | -0.040 | -0.056 | -0.011
  d5 ZT0-6 +26.3 | +0.032 | +0.028 | +0.041     d6 ZT0-6 +29.9 | -0.038 | +0.028 | +0.043
  d7 ZT0-6 +33.6 | +0.002 | +0.000 | +0.013     d8 ZT0-6 +57.6 | -0.104 | -0.111 | +0.019
POST HOC (descriptive, not in the PRESPEC): sleep placebo day 2 - days 0-1 median +28.0 min, p 0.0073 |
  day-4 increase minus placebo median +11.6 min, paired p 0.094
```
VERDICTS by the pre-specified rules:
- Stage 0 PASS (+42 min).
- All three primary readouts NOT SUPPORTED. Wording: "No change of c_obs, cycle flux or count EPR beyond day-to-day
  variation was detected in the first 6 h of recovery (median Delta c_obs -0.062 [-0.151, +0.004], J -0.056 cycles/h
  [-0.194, +0.056], EPR +0.037 nats/h [-0.003, +0.152])." EPR's nominal p 0.024 does not survive Holm (0.071), and its
  placebo change (+0.091) is LARGER than the recovery change; its time-trend check is borderline (0.059).
- Transitions: by the pre-specified Holm, "sleep pressure changes" W->S (x1.49), W->Q (x1.17) and S->W (x0.75).
  Recovery flies fall asleep from active wake more readily, enter quiet wake more readily, and wake into active wake
  less readily.

CRITICAL CAVEAT (found in the pre-specified descriptive rows, confirmed by one post hoc line):
- ZT0-6 sleep drifts upward through the recording: +26 to +58 min on days 5-8, and day 8 (+58) exceeds day 4 (+42).
- The same placebo logic the primary readouts use, applied post hoc to sleep, is itself significant (day 2 vs days
  0-1: +28 min, p 0.007). The day-4 increase beyond it is +11.6 min (p 0.094).
- So the Stage 0 pass CANNOT be read as a deprivation rebound distinct from the time trend. The gate as specified had
  no placebo; that was a design gap (mine, in writing the PRESPEC to the brief).
- The transition ratios also had no placebo control. They may reflect the same upward drift in sleep propensity
  rather than deprivation. They must NOT be reported as deprivation effects without that caveat.
- Reportable reading: in these 26 flies, a 9-h night-time deprivation produced no change in the cycle's direction,
  flux or dissipation beyond day-to-day variation; even the sleep rebound itself is not separable from a steady upward
  trend in morning sleep across the 9.5-day recording.

## 20x post hoc P1-P2 (2026-11-01; POST HOC, specified after the 20x results; note 3e6f979, code 78748d8; planB_20x_posthoc_P1P2.log, planB_20x_posthoc.json)
Narrowest block: one fly's ZT0-6 on day 4. Placebo arm = day 2 vs days 0-1, same flies, same window and estimator.
```
P1 transition | recovery/baseline [CI] | placebo day2/days0-1 [CI] | log contrast [CI] p Holm
  Q->W | 0.94 [0.83, 1.05] | 1.02 [0.92, 1.14] | -0.09 [-0.22, +0.06] p 0.217 Holm 0.868
  Q->S | 1.16 [0.98, 1.34] | 0.97 [0.78, 1.18] | +0.18 [+0.01, +0.36] p 0.033 Holm 0.198
  W->Q | 1.17 [1.05, 1.31] | 1.20 [1.04, 1.37] | -0.02 [-0.19, +0.17] p 0.857 Holm 1.000
  W->S | 1.49 [1.24, 1.79] | 1.26 [1.00, 1.54] | +0.17 [-0.06, +0.45] p 0.160 Holm 0.800
  S->Q | 0.90 [0.71, 1.12] | 0.90 [0.72, 1.15] | -0.00 [-0.28, +0.28] p 0.962 Holm 1.000
  S->W | 0.75 [0.58, 0.95] | 0.86 [0.68, 1.07] | -0.13 [-0.39, +0.15] p 0.346 Holm 1.000
P2 median active fraction of non-S minutes
  ZT0-6   baseline 0.067 | Delta -0.0167 p 0.0048 | Delta_pl -0.0167 | Delta vs Delta_pl paired p 0.196 -> no change beyond day-to-day variation
  ZT6-12  baseline 0.050 | Delta +0.0000 p 0.42   | Delta_pl +0.0000 | paired p 0.088 (descriptive)
  ZT12-24 baseline 0.200 | Delta -0.0500 p 8.5e-05 | Delta_pl -0.0333 | paired p 0.033 (descriptive)
```
- P1: NONE of the six transitions changes beyond the day-to-day trend. All read "not separable from the trend in morning
  sleep". The 20x transition findings do not survive a placebo arm:
  - W->Q changes as much between baseline days (x1.20) as after deprivation (x1.17);
  - W->S (x1.49 vs x1.26) and S->W (x0.75 vs x0.86) move in the same direction on the placebo days;
  - Q->S has the only nominal contrast (+0.18, p 0.033, CI just excluding 0), which does not survive Holm (0.198).
- P2: wake activity at ZT0-6 is lower after deprivation (-0.017, p 0.005), but by exactly as much between baseline days
  (placebo -0.017): no change beyond day-to-day variation. Descriptively, day 4 ZT12-24 drops more than its placebo
  (-0.050 vs -0.033, paired p 0.033); one descriptive line among several, not a result.
- Combined reading of 20x: in this dataset, no measured property of the observed cycle (direction, flux, dissipation,
  transition probabilities, wake activity level) changes after 9 h of deprivation beyond the drift seen between
  baseline days. A parallel undeprived control cohort is needed to test the homeostat on the cycle.

## PRESPEC 20u — reversible-null control of the per0 DD wake-hazard rise: VERDICT NULL INADEQUATE (2026-11-01; PRESPEC ccae0e8 + amendment A 3f3f2cd; code cad83a4 + checkpointing aa026cd; planB_20u_run_part1.log + planB_20u_run.log, planB_20u_run.json)
Stage 1 ran in two parts: the first background job was stopped by the harness's 2-h limit after 23 of 25 flies were
checkpointed, and the resumed job fitted the last 2. Stage 1 total ~2.5 h, inside the 6-h box. The earlier crash
(underflow) is noted above.
```
selected orders (Ka, Ki): (3,3) in 15 flies, (2,3) in 10; runs per fly 4,648-37,895; EM iterations 82-1999
laziness: smallest eigenvalue of the selected P <= 0 in 24 of 25 flies (range -0.296 to -0.010; fly 16 +0.054)
stage 1b adequacy: wake-bout count ratio median 1.225 [IQR 1.136, 1.278] | total sleep ratio median 1.015 [IQR 0.999, 1.031]
                   | death-excluded sims 0/500 -> FAIL
VERDICT: NULL INADEQUATE (stage 2 not run, per PRESPEC 20u)
```
- By the pre-specified rule the test says NOTHING either way about the per0 DD hazard rise: the best-fitting reversible
  aggregated HMM (K <= 6) does not reproduce the bout process. It makes 22.5% MORE wake bouts than observed (20r-B's
  renewal null: ~15% more), while total sleep is right (x1.015). Real flies' wake is broken into sleep less often than
  reversible dynamics with the same 1-s statistics predict.
- That over-fragmentation is itself an observation, but NOT a test of reversibility. 15/25 flies selected the
  largest allowed order (3+3), so the model class was at its ceiling; a larger reversible model might fit. Kmax was
  fixed in advance and is not raised (no follow-up without a new PRESPEC).
- Laziness failed in 24/25 flies: the fitted 1-s chains have negative eigenvalues (fast flip-flop between hidden
  states). The discrete-time complete-monotonicity premise of the motivating theorem therefore does not hold for these
  fits either. It was motivation only, not the test.
- C8 stays as written ("supported but bin-sensitive"; 20r-B inconclusive). 20u adds no evidence for or against.
