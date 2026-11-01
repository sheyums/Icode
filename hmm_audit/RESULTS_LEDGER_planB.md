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
