# Results ledger (append-only). Code commit: fd9dcc9. PRESPEC commit: 37d2684. Date 2026-10-19.

## A (planned) emission degeneracy, L1 = sum|B1-B2|
- CS 594 windows/25 flies: L1<1e-6 0.678; <0.05 0.690; <0.2 0.690. Median B = [0 1; 0 1; 1 0].
- fumin 1466 windows/53 flies: L1<1e-6 0.993; <0.05 0.993; <0.2 0.993.
- Unplanned look at the rest: H1 emits like H3 (inactive) instead: CS 163/191 rounded, 31 exact; fumin 9/10 rounded, 8 exact.

## B (planned) identified set of Xn_pred / Yn_pred, exactly H1==H2 windows (inner bounds, grid)
- CS n=403: Xn undetermined 0.993, Yn 0.998; median interval Xn [-0.996,+0.993], Yn [-1,+1].
- fumin n=1456: Xn undetermined 1.000, Yn 1.000; median interval Xn [-0.997,+0.997], Yn [-1,+1].
- Spot checks: logL equal to 1e-6 under P and extreme P' (6/6 windows).
- Unplanned: fixed uniform pi0, extreme P' vs P logL difference in [-0.031, +0.683] nats (CS), 0.000 (fumin).
- Unplanned B13 (exact H1==H3 windows): CS n=31 Xn/Yn undetermined 1.000/1.000; fumin n=8 1.000/0.875.

## C (planned) Viterbi vs FFBS mean, within-window, raw R (no smoothing), M=200, 0 FFBS failures
- CS: Viterbi N32 1957 N312 1096 N23 1333 N213 1755 (Xn -0.282, Yn +0.137); FFBS 1623.0 1438.5 1102.6 1994.5 (Xn -0.060, Yn +0.288).
- fumin: Viterbi 5239 155 854 3433 (Xn -0.943, Yn +0.602); FFBS 4085.5 1312.4 399.7 3892.4 (Xn -0.514, Yn +0.814).
- Per-window Viterbi-FFBS N312: CS mean -0.58 median 0; fumin mean -0.79 median -0.64.

## 2026-10-19 — PRESPEC 2 (468b1bb); code committed with this entry
### A/B sss, per0
    
    === sss: 593 accepted windows, 57 flies
      A: L1 < 1e-06: 0.988  (586)
      A: L1 < 0.05: 0.988  (586)
      A: L1 < 0.2: 0.988  (586)
      A: median B rows: [[0.0, 1.0], [0.0, 1.0], [1.0, 0.0]]
      B: exactly-degenerate windows analysed: 586
      B: fitted Xn_pred  median -0.470;  identified interval median [-0.997, +0.997]
      B: fitted Yn_pred  median +1.000;  identified interval median [-1.000, +1.000]
      B: Xn undefined in all models: 0.000; Yn: 0.000
      B: Xn direction UNDETERMINED (interval contains 0): 1.000
      B: Yn direction UNDETERMINED (interval contains 0): 1.000
      B: Xn interval width median 1.995; Yn 2.000
      spot fly1 win2: logL P -26.896860 | P' -26.896860 | Xn_pred -0.470 -> +0.997
      spot fly2 win1: logL P -43.094730 | P' -43.094730 | Xn_pred -0.471 -> +0.997
      spot fly3 win1: logL P -8.952032 | P' -8.952032 | Xn_pred -0.406 -> +0.996
    
    === per0: 544 accepted windows, 25 flies
      A: L1 < 1e-06: 0.570  (310)
      A: L1 < 0.05: 0.599  (326)
      A: L1 < 0.2: 0.601  (327)
      A: median B rows: [[0.0, 1.0], [0.0, 1.0], [1.0, 0.0]]
      B: exactly-degenerate windows analysed: 310
      B: fitted Xn_pred  median +0.035;  identified interval median [-0.997, +0.994]
      B: fitted Yn_pred  median +0.985;  identified interval median [-1.000, +1.000]
      B: Xn undefined in all models: 0.000; Yn: 0.000
      B: Xn direction UNDETERMINED (interval contains 0): 1.000
      B: Yn direction UNDETERMINED (interval contains 0): 0.997
      B: Xn interval width median 1.989; Yn 2.000
      spot fly3 win1: logL P -125.517071 | P' -125.517071 | Xn_pred -1.000 -> +0.833
      spot fly6 win1: logL P -105.063022 | P' -105.063022 | Xn_pred -1.000 -> +1.000
      spot fly7 win1: logL P -119.261848 | P' -119.261848 | Xn_pred -0.263 -> +0.949
### C sss, per0 (0 FFBS failures)
    sss: Viterbi N32 1559 N312 64 N23 191 N213 971 Xn -0.921 Yn +0.671 | FFBS 1204.6 422.2 94.3 1071.5 Xn -0.481 Yn +0.838
    per0: Viterbi N32 1865 N312 1668 N23 2198 N213 1525 Xn -0.056 Yn -0.181 | FFBS 1676.2 1864.9 1896.3 1834.8 Xn +0.053 Yn -0.016
### D all genotypes
    
    CS     H1==H2: 403 windows
       Xn           invariant 0.000  sign determined 0.000  median width 1.99  undefined 3  | window-mean identified set [-0.992, +0.985]
       Yn           invariant 0.000  sign determined 0.002  median width 2.00  undefined 0  | window-mean identified set [-0.997, +0.990]
       Delta=Xn-Yn  invariant 0.000  sign determined 0.000  median width 3.77  undefined 3  | window-mean identified set [-1.882, +1.832]
       Sigma=Xn+Yn  invariant 0.058  sign determined 0.000  median width 3.95  undefined 3  | window-mean identified set [-1.732, +1.715]
       spot fly1 win1: fit: Xn -0.04 Yn +1.00 D -1.04 logL -113.6194 | maxD: Xn +0.96 Yn -0.90 D +1.86 logL -113.6194 | minD: Xn -0.97 Yn +0.90 D -1.87 logL -113.6194
       spot fly2 win1: fit: Xn -0.87 Yn +1.00 D -1.87 logL -50.1711 | maxD: Xn +0.98 Yn -1.00 D +1.98 logL -50.1711 | minD: Xn -0.97 Yn +1.00 D -1.97 logL -50.1711
       spot fly3 win1: fit: Xn -1.00 Yn +1.00 D -2.00 logL -95.2240 | maxD: Xn +1.00 Yn -1.00 D +2.00 logL -95.2240 | minD: Xn -1.00 Yn +1.00 D -2.00 logL -95.2240
    
    CS     H1==H3: 31 windows
       Xn           invariant 0.000  sign determined 0.000  median width 2.00  undefined 0  | window-mean identified set [-0.999, +0.999]
       Yn           invariant 0.000  sign determined 0.000  median width 1.99  undefined 0  | window-mean identified set [-0.995, +0.995]
       Delta=Xn-Yn  invariant 0.000  sign determined 0.000  median width 3.91  undefined 0  | window-mean identified set [-1.953, +1.953]
       Sigma=Xn+Yn  invariant 0.097  sign determined 0.000  median width 3.20  undefined 0  | window-mean identified set [-1.299, +1.286]
    
    fumin  H1==H2: 1456 windows
       Xn           invariant 0.000  sign determined 0.000  median width 1.99  undefined 0  | window-mean identified set [-0.997, +0.996]
       Yn           invariant 0.000  sign determined 0.000  median width 2.00  undefined 0  | window-mean identified set [-1.000, +1.000]
       Delta=Xn-Yn  invariant 0.000  sign determined 0.000  median width 3.78  undefined 0  | window-mean identified set [-1.886, +1.885]
       Sigma=Xn+Yn  invariant 0.021  sign determined 0.000  median width 2.72  undefined 0  | window-mean identified set [-1.332, +1.316]
       spot fly1 win5: fit: Xn -0.55 Yn +1.00 D -1.55 logL -40.8258 | maxD: Xn +0.99 Yn -0.93 D +1.92 logL -40.8258 | minD: Xn -0.99 Yn +0.96 D -1.95 logL -40.8258
       spot fly2 win12: fit: Xn -0.52 Yn +nan D +nan logL -9.0906 | maxD: Xn +0.88 Yn +nan D +nan logL -9.0906 | minD: Xn -0.88 Yn +nan D +nan logL -9.0906
       spot fly3 win2: fit: Xn -0.42 Yn +1.00 D -1.42 logL -26.5060 | maxD: Xn +0.83 Yn -1.00 D +1.83 logL -26.5060 | minD: Xn -0.87 Yn +1.00 D -1.87 logL -26.5060
    
    fumin  H1==H3: 8 windows
       Xn           invariant 0.000  sign determined 0.000  median width 2.00  undefined 0  | window-mean identified set [-1.000, +1.000]
       Yn           invariant 0.000  sign determined 0.000  median width 1.99  undefined 1  | window-mean identified set [-0.997, +0.997]
       Delta=Xn-Yn  invariant 0.000  sign determined 0.000  median width 3.94  undefined 1  | window-mean identified set [-1.973, +1.976]
       Sigma=Xn+Yn  invariant 0.286  sign determined 0.000  median width 1.63  undefined 1  | window-mean identified set [-0.703, +0.667]
    
    sss    H1==H2: 586 windows
       Xn           invariant 0.000  sign determined 0.000  median width 1.99  undefined 0  | window-mean identified set [-0.997, +0.997]
       Yn           invariant 0.000  sign determined 0.000  median width 2.00  undefined 0  | window-mean identified set [-1.000, +1.000]
       Delta=Xn-Yn  invariant 0.000  sign determined 0.000  median width 3.77  undefined 0  | window-mean identified set [-1.890, +1.889]
       Sigma=Xn+Yn  invariant 0.036  sign determined 0.000  median width 2.84  undefined 0  | window-mean identified set [-1.376, +1.363]
       spot fly1 win2: fit: Xn -0.42 Yn +1.00 D -1.42 logL -26.8969 | maxD: Xn +0.88 Yn -1.00 D +1.88 logL -26.8969 | minD: Xn -0.89 Yn +1.00 D -1.89 logL -26.8969
       spot fly2 win1: fit: Xn -0.49 Yn +1.00 D -1.49 logL -43.0947 | maxD: Xn +0.93 Yn -0.97 D +1.89 logL -43.0947 | minD: Xn -0.94 Yn +0.94 D -1.88 logL -43.0947
       spot fly3 win1: fit: Xn -0.43 Yn +nan D +nan logL -8.9520 | maxD: Xn +0.84 Yn +nan D +nan logL -8.9520 | minD: Xn -0.85 Yn +nan D +nan logL -8.9520
    
    sss    H1==H3: 7 windows
       Xn           invariant 0.000  sign determined 0.000  median width 2.00  undefined 0  | window-mean identified set [-1.000, +1.000]
       Yn           invariant 0.000  sign determined 0.000  median width 2.00  undefined 1  | window-mean identified set [-0.998, +0.998]
       Delta=Xn-Yn  invariant 0.000  sign determined 0.000  median width 3.95  undefined 1  | window-mean identified set [-1.957, +1.958]
       Sigma=Xn+Yn  invariant 0.000  sign determined 0.000  median width 2.37  undefined 1  | window-mean identified set [-1.229, +1.194]
    
    per0   H1==H2: 310 windows
       Xn           invariant 0.000  sign determined 0.000  median width 1.99  undefined 0  | window-mean identified set [-0.993, +0.986]
       Yn           invariant 0.000  sign determined 0.003  median width 2.00  undefined 0  | window-mean identified set [-0.993, +0.978]
       Delta=Xn-Yn  invariant 0.000  sign determined 0.000  median width 3.80  undefined 0  | window-mean identified set [-1.900, +1.845]
       Sigma=Xn+Yn  invariant 0.110  sign determined 0.013  median width 3.93  undefined 0  | window-mean identified set [-1.651, +1.616]
       spot fly1 win2: fit: Xn +0.07 Yn -0.80 D +0.87 logL -68.5570 | maxD: Xn +0.99 Yn -0.98 D +1.97 logL -68.5570 | minD: Xn -0.99 Yn +0.98 D -1.97 logL -68.5570
       spot fly2 win2: fit: Xn +0.44 Yn +0.57 D -0.13 logL -79.0173 | maxD: Xn +0.95 Yn -0.99 D +1.94 logL -79.0173 | minD: Xn -0.97 Yn +0.99 D -1.96 logL -79.0173
       spot fly3 win1: fit: Xn -1.00 Yn +1.00 D -2.00 logL -125.5171 | maxD: Xn +1.00 Yn -1.00 D +2.00 logL -125.5171 | minD: Xn -1.00 Yn +1.00 D -2.00 logL -125.5171
    
    per0   H1==H3: 75 windows
       Xn           invariant 0.000  sign determined 0.000  median width 2.00  undefined 0  | window-mean identified set [-0.999, +0.997]
       Yn           invariant 0.000  sign determined 0.000  median width 2.00  undefined 0  | window-mean identified set [-0.996, +0.995]
       Delta=Xn-Yn  invariant 0.000  sign determined 0.000  median width 3.99  undefined 0  | window-mean identified set [-1.976, +1.974]
       Sigma=Xn+Yn  invariant 0.467  sign determined 0.000  median width 0.08  undefined 0  | window-mean identified set [-0.647, +0.641]

## 2026-10-19 — PRESPEC 3 (7e19368), analysis E: ABANDONED before any result
The script crashed on a shape bug before printing any E result (log showed only a
re-run of D, identical to the D entry above). The user then clarified that the
comparison of interest is N312 vs N32 and N213 vs N23 (= Xn, Yn, already in B/D),
not N312 vs N213. E was not re-run. No E result exists.

## 2026-10-19 — PRESPEC 4 (9c8c1ac): per-fly scatter (plot_scatter.py, committed with this entry)
    CS     flies  25 | FFBS Xn median +0.02 Yn +0.68 | Viterbi Xn -0.29 Yn +0.63 | median range Xn [-0.99,+0.98] Yn [-1.00,+1.00]
    per0   flies  25 | FFBS Xn median -0.03 Yn +0.40 | Viterbi Xn -0.11 Yn +0.27 | median range Xn [-1.00,+0.99] Yn [-0.99,+0.99]
    fumin  flies  53 | FFBS Xn median -0.52 Yn +0.87 | Viterbi Xn -1.00 Yn +0.68 | median range Xn [-1.00,+1.00] Yn [-1.00,+1.00]
    sss    flies  54 | FFBS Xn median -0.50 Yn +0.95 | Viterbi Xn -1.00 Yn +0.87 | median range Xn [-1.00,+1.00] Yn [-1.00,+1.00]
    sss: 3 of 57 flies dropped (no exact H1==H2 window with defined Xn and Yn), 575 of 586 windows plotted.
    CORRECTION to the line above: the 3 sss flies (5, 6, 51; 11 windows) have qualifying windows but zero posterior wake->sleep routes (N23+N213 = 0), so Yn is undefined.

## 2026-10-19 — descriptive (unplanned, no PRESPEC): runs per accepted window in raw R, and model index counts in C
    CS     runs (bouts) per window: median 11, IQR 8-14, min 2; frames median 14212 | selected m: {2: 2, 3: 730, 4: 6}
    per0   runs (bouts) per window: median 15, IQR 11-18, min 2; frames median 14352 | selected m: {3: 628, 4: 2}
    fumin  runs (bouts) per window: median 6, IQR 4-10, min 2; frames median 13669 | selected m: {2: 298, 3: 1466, 4: 5}
    sss    runs (bouts) per window: median 4, IQR 2-8, min 1; frames median 12379 | selected m: {2: 241, 3: 828, 4: 3, 5: 1}

## 2026-10-19 — PRESPEC 5 (eb4076b), pooled.py (committed 0822229) — PRIMARY, whole record
    CS: flies 25, sleep bouts 3410, wake bouts 3413
    per0: flies 25, sleep bouts 3986, wake bouts 3988
    fumin: flies 53, sleep bouts 5263, wake bouts 5269
    sss: flies 55, sleep bouts 1486, wake bouts 1522
    
    === P1: phases per bout class (discrete, left-truncated hyperexponential) ===
    
    CS sleep n=3410  CV-selected K=3  (AICc-best K=3, BIC-best K=3)
       K=1 logL   -27430.96  CV logL   -27443.12  AICc    54863.9  BIC    54870.1 | q [1.0] means(s) [1446.0] min branch 3410.0
       K=2 logL   -27275.16  CV logL   -27288.15  AICc    54556.3  BIC    54574.7 | q [0.19, 0.81] means(s) [437.0, 1682.0] min branch 647.3
       K=3 logL   -27250.32  CV logL   -27263.13  AICc    54510.7  BIC    54541.3 | q [0.137, 0.752, 0.111] means(s) [398.0, 1390.0, 3117.0] min branch 378.1
       K=4 logL   -27250.31  CV logL   -27265.16  AICc    54514.7  BIC    54557.6 | q [0.0, 0.137, 0.752, 0.111] means(s) [337.0, 398.0, 1390.0, 3117.0] min branch 0.6
    
    CS wake  n=3413  CV-selected K=4  (AICc-best K=3, BIC-best K=3)
       K=1 logL   -26952.86  CV logL   -26976.26  AICc    53907.7  BIC    53913.9 | q [1.0] means(s) [991.0] min branch 3413.0
       K=2 logL   -26631.83  CV logL   -26651.07  AICc    53269.7  BIC    53288.1 | q [0.789, 0.211] means(s) [549.0, 2639.0] min branch 721.6
       K=3 logL   -26476.20  CV logL   -26496.37  AICc    52962.4  BIC    52993.1 | q [0.048, 0.771, 0.181] means(s) [5.0, 620.0, 2837.0] min branch 164.1
       K=4 logL   -26474.30  CV logL   -26495.53  AICc    52962.6  BIC    53005.5 | q [0.042, 0.018, 0.765, 0.175] means(s) [5.0, 36.0, 636.0, 2878.0] min branch 61.5
    
    per0 sleep n=3986  CV-selected K=2  (AICc-best K=3, BIC-best K=3)
       K=1 logL   -31656.63  CV logL   -31663.53  AICc    63315.3  BIC    63321.6 | q [1.0] means(s) [1334.0] min branch 3986.0
       K=2 logL   -31618.22  CV logL   -31632.13  AICc    63242.5  BIC    63261.3 | q [0.056, 0.944] means(s) [360.0, 1392.0] min branch 223.8
       K=3 logL   -31609.84  CV logL   -31641.13  AICc    63229.7  BIC    63261.1 | q [0.05, 0.93, 0.021] means(s) [354.0, 1340.0, 3436.0] min branch 82.7
       K=4 logL   -31609.84  CV logL   -31637.40  AICc    63233.7  BIC    63277.7 | q [0.0, 0.05, 0.93, 0.021] means(s) [302.0, 354.0, 1340.0, 3436.0] min branch 0.0
    
    per0 wake  n=3988  CV-selected K=2  (AICc-best K=2, BIC-best K=2)
       K=1 logL   -29561.77  CV logL   -29573.52  AICc    59125.5  BIC    59131.8 | q [1.0] means(s) [611.0] min branch 3988.0
       K=2 logL   -29317.73  CV logL   -29336.90  AICc    58641.5  BIC    58660.3 | q [0.079, 0.921] means(s) [10.0, 663.0] min branch 316.3
       K=3 logL   -29317.40  CV logL   -29342.38  AICc    58644.8  BIC    58676.3 | q [0.079, 0.919, 0.002] means(s) [10.0, 659.0, 2434.0] min branch 7.0
       K=4 logL   -29317.34  CV logL   -29342.10  AICc    58648.7  BIC    58692.7 | q [0.047, 0.033, 0.918, 0.002] means(s) [8.0, 13.0, 660.0, 2454.0] min branch 6.8
    
    fumin sleep n=5263  CV-selected K=4  (AICc-best K=3, BIC-best K=2)
       K=1 logL   -33479.44  CV logL   -33502.47  AICc    66960.9  BIC    66967.4 | q [1.0] means(s) [512.0] min branch 5263.0
       K=2 logL   -33431.02  CV logL   -33454.41  AICc    66868.0  BIC    66887.7 | q [0.727, 0.273] means(s) [459.0, 653.0] min branch 1438.4
       K=3 logL   -33428.33  CV logL   -33454.54  AICc    66866.7  BIC    66899.5 | q [0.02, 0.788, 0.192] means(s) [319.0, 475.0, 687.0] min branch 106.6
       K=4 logL   -33428.33  CV logL   -33454.20  AICc    66870.7  BIC    66916.6 | q [0.021, 0.784, 0.098, 0.097] means(s) [320.0, 475.0, 661.0, 708.0] min branch 111.3
    
    fumin wake  n=5269  CV-selected K=3  (AICc-best K=3, BIC-best K=3)
       K=1 logL   -46998.30  CV logL   -47001.21  AICc    93998.6  BIC    94005.2 | q [1.0] means(s) [2753.0] min branch 5269.0
       K=2 logL   -46760.80  CV logL   -46764.82  AICc    93527.6  BIC    93547.3 | q [0.676, 0.324] means(s) [1514.0, 5335.0] min branch 1707.7
       K=3 logL   -46749.07  CV logL   -46756.02  AICc    93508.1  BIC    93541.0 | q [0.006, 0.675, 0.319] means(s) [7.0, 1544.0, 5366.0] min branch 32.3
       K=4 logL   -46748.89  CV logL   -46756.33  AICc    93511.8  BIC    93557.8 | q [0.0, 0.006, 0.675, 0.319] means(s) [2.0, 8.0, 1544.0, 5366.0] min branch 2.0
    
    sss sleep n=1486  CV-selected K=2  (AICc-best K=2, BIC-best K=2)
       K=1 logL    -9641.30  CV logL    -9660.92  AICc    19284.6  BIC    19289.9 | q [1.0] means(s) [541.0] min branch 1486.0
       K=2 logL    -9559.98  CV logL    -9574.58  AICc    19126.0  BIC    19141.9 | q [0.681, 0.319] means(s) [429.0, 780.0] min branch 474.7
       K=3 logL    -9558.90  CV logL    -9574.94  AICc    19127.8  BIC    19154.3 | q [0.307, 0.488, 0.205] means(s) [388.0, 503.0, 861.0] min branch 305.1
       K=4 logL    -9558.90  CV logL    -9574.59  AICc    19131.9  BIC    19168.9 | q [0.0, 0.307, 0.488, 0.205] means(s) [302.0, 388.0, 503.0, 861.0] min branch 0.0
    
    sss wake  n=1522  CV-selected K=4  (AICc-best K=3, BIC-best K=3)
       K=1 logL   -13933.37  CV logL   -13942.94  AICc    27868.7  BIC    27874.1 | q [1.0] means(s) [3481.0] min branch 1522.0
       K=2 logL   -13766.10  CV logL   -13776.80  AICc    27538.2  BIC    27554.2 | q [0.418, 0.582] means(s) [756.0, 5441.0] min branch 636.6
       K=3 logL   -13735.35  CV logL   -13748.90  AICc    27480.7  BIC    27507.3 | q [0.037, 0.422, 0.541] means(s) [17.0, 968.0, 5677.0] min branch 56.3
       K=4 logL   -13734.10  CV logL   -13748.02  AICc    27482.3  BIC    27519.5 | q [0.022, 0.027, 0.419, 0.532] means(s) [11.0, 63.0, 1030.0, 5730.0] min branch 33.7
    
    === P2: consecutive-bout coupling (within-window centred ranks) ===
    CS     pairs SW 3255 WS 3276 | rho_SW +0.034 CI [-0.007  0.072] perm p 0.0550 | rho_WS -0.020 CI [-0.06   0.022] perm p 0.3150 | A +0.054 CI [-0.008  0.116]
    per0   pairs SW 3934 WS 3945 | rho_SW +0.086 CI [0.039 0.132] perm p 0.0000 | rho_WS +0.094 CI [0.044 0.139] perm p 0.0000 | A -0.008 CI [-0.068  0.057]
    fumin  pairs SW 4276 WS 4228 | rho_SW +0.026 CI [-0.012  0.064] perm p 0.1255 | rho_WS +0.013 CI [-0.013  0.038] perm p 0.4620 | A +0.013 CI [-0.036  0.065]
    sss    pairs SW 1035 WS 990 | rho_SW +0.004 CI [-0.055  0.067] perm p 0.9150 | rho_WS -0.025 CI [-0.093  0.04 ] perm p 0.4765 | A +0.029 CI [-0.061  0.124]
    
    Genotype differences in A (bootstrap CI of difference):
       A[CS] - A[per0]: +0.061 CI [-0.025  0.148]
       A[CS] - A[fumin]: +0.038 CI [-0.04   0.118]
       A[CS] - A[sss]: +0.024 CI [-0.09   0.136]
       A[per0] - A[fumin]: -0.021 CI [-0.096  0.056]
       A[per0] - A[sss]: -0.039 CI [-0.151  0.076]
       A[fumin] - A[sss]: -0.017 CI [-0.123  0.089]

## 2026-10-19 — PRESPEC 5b run 1 CRASHED (seed bug: negative fold code as RNG seed). Only output: crossing fractions (kept in pooled_LD_crash1.log). Fix: seed offset +20 in an LD-only job function; rerun.
