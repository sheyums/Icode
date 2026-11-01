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
