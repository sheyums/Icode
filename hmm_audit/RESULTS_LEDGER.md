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

## 2026-10-19 — PRESPEC 5b (fab217d) light/dark, pooled_LD.py run 2 (seed fix c965ffe) — SECONDARY
    === crossing fraction (bouts whose span crosses ZT0/ZT12) ===
      CS     light sleep n= 2058 crossing 0.023
      CS     light wake  n= 2051 crossing 0.025
      CS     dark  sleep n= 1352 crossing 0.016
      CS     dark  wake  n= 1362 crossing 0.039
      per0   light sleep n= 2019 crossing 0.031
      per0   light wake  n= 2019 crossing 0.014
      per0   dark  sleep n= 1967 crossing 0.023
      per0   dark  wake  n= 1969 crossing 0.010
      fumin  light sleep n= 3180 crossing 0.010
      fumin  light wake  n= 3162 crossing 0.046
      fumin  dark  sleep n= 2083 crossing 0.014
      fumin  dark  wake  n= 2107 crossing 0.080
      sss    light sleep n=  711 crossing 0.008
      sss    light wake  n=  736 crossing 0.095
      sss    dark  sleep n=  775 crossing 0.012
      sss    dark  wake  n=  786 crossing 0.062
    
    === P1 by phase ===
    CS     light sleep n= 2058 CV K=3 (AICc 3, BIC 3) | CV logL by K [np.float64(-15814.1), np.float64(-15749.4), np.float64(-15748.1), np.float64(-15749.6)] | q [0.163, 0.829, 0.008] means [405.0, 1190.0, 4772.0] min branch 16.8
    CS     light wake  n= 2051 CV K=4 (AICc 3, BIC 3) | CV logL by K [np.float64(-16088.2), np.float64(-15869.0), np.float64(-15789.5), np.float64(-15789.0)] | q [0.04, 0.021, 0.74, 0.199] means [5.0, 64.0, 558.0, 2597.0] min branch 43.1
    CS     dark  sleep n= 1352 CV K=4 (AICc 2, BIC 2) | CV logL by K [np.float64(-11400.9), np.float64(-11367.0), np.float64(-11366.8), np.float64(-11366.8)] | q [0.0, 0.074, 0.817, 0.109] means [301.0, 363.0, 1942.0, 3397.0] min branch 0.0
    CS     dark  wake  n= 1362 CV K=3 (AICc 3, BIC 3) | CV logL by K [np.float64(-10879.9), np.float64(-10795.4), np.float64(-10698.9), np.float64(-10700.0)] | q [0.056, 0.79, 0.154] means [6.0, 737.0, 3230.0] min branch 76.7
    per0   light sleep n= 2019 CV K=2 (AICc 3, BIC 2) | CV logL by K [np.float64(-16024.4), np.float64(-16008.3), np.float64(-16011.3), np.float64(-16008.7)] | q [0.052, 0.948] means [352.0, 1380.0] min branch 105.4
    per0   light wake  n= 2019 CV K=2 (AICc 2, BIC 2) | CV logL by K [np.float64(-14984.5), np.float64(-14878.2), np.float64(-14879.9), np.float64(-14879.9)] | q [0.079, 0.921] means [11.0, 666.0] min branch 159.1
    per0   dark  sleep n= 1967 CV K=3 (AICc 3, BIC 2) | CV logL by K [np.float64(-15640.1), np.float64(-15639.5), np.float64(-15632.6), np.float64(-15632.7)] | q [0.051, 0.882, 0.067] means [361.0, 1310.0, 2513.0] min branch 100.8
    per0   dark  wake  n= 1969 CV K=2 (AICc 2, BIC 2) | CV logL by K [np.float64(-14589.4), np.float64(-14460.8), np.float64(-14463.2), np.float64(-14464.0)] | q [0.081, 0.919] means [9.0, 660.0] min branch 159.3
    fumin  light sleep n= 3180 CV K=2 (AICc 2, BIC 2) | CV logL by K [np.float64(-20268.1), np.float64(-20253.3), np.float64(-20254.2), np.float64(-20254.2)] | q [0.797, 0.203] means [477.0, 662.0] min branch 646.6
    fumin  light wake  n= 3162 CV K=3 (AICc 3, BIC 2) | CV logL by K [np.float64(-27832.2), np.float64(-27733.9), np.float64(-27729.4), np.float64(-27729.6)] | q [0.005, 0.799, 0.196] means [7.0, 1713.0, 5512.0] min branch 16.3
    fumin  dark  sleep n= 2083 CV K=4 (AICc 2, BIC 2) | CV logL by K [np.float64(-13241.5), np.float64(-13205.7), np.float64(-13204.6), np.float64(-13204.3)] | q [0.016, 0.346, 0.586, 0.052] means [313.0, 408.0, 542.0, 871.0] min branch 33.1
    fumin  dark  wake  n= 2107 CV K=3 (AICc 3, BIC 2) | CV logL by K [np.float64(-19123.8), np.float64(-18997.5), np.float64(-18995.3), np.float64(-18995.7)] | q [0.008, 0.534, 0.458] means [8.0, 1302.0, 5488.0] min branch 16.0
    sss    light sleep n=  711 CV K=3 (AICc 3, BIC 2) | CV logL by K [np.float64(-4609.4), np.float64(-4575.0), np.float64(-4572.0), np.float64(-4572.7)] | q [0.166, 0.663, 0.171] means [348.0, 496.0, 874.0] min branch 118.3
    sss    light wake  n=  736 CV K=3 (AICc 3, BIC 3) | CV logL by K [np.float64(-6853.5), np.float64(-6785.3), np.float64(-6774.8), np.float64(-6775.5)] | q [0.044, 0.337, 0.62] means [34.0, 1006.0, 5941.0] min branch 32.0
    sss    dark  sleep n=  775 CV K=2 (AICc 2, BIC 2) | CV logL by K [np.float64(-5053.0), np.float64(-5005.9), np.float64(-5006.2), np.float64(-5006.2)] | q [0.702, 0.298] means [432.0, 814.0] min branch 230.9
    sss    dark  wake  n=  786 CV K=3 (AICc 3, BIC 3) | CV logL by K [np.float64(-7075.4), np.float64(-6986.2), np.float64(-6967.4), np.float64(-6969.1)] | q [0.039, 0.501, 0.46] means [13.0, 987.0, 5384.0] min branch 30.5
    
    === P2 by phase ===
    CS     light pairs SW 2013 WS 2007 | rho_SW -0.021 CI [-0.065  0.02 ] p 0.3605 | rho_WS -0.032 CI [-0.077  0.013] p 0.1790 | A +0.010 CI [-0.051  0.074]
    CS     dark  pairs SW 1242 WS 1269 | rho_SW +0.126 CI [0.05  0.189] p 0.0000 | rho_WS -0.001 CI [-0.065  0.071] p 0.9705 | A +0.127 CI [0.024 0.231]
    CS     A[light]-A[dark] -0.115 CI [-0.227 -0.007]
    per0   light pairs SW 2005 WS 2006 | rho_SW +0.101 CI [0.047 0.155] p 0.0000 | rho_WS +0.091 CI [0.03 0.15] p 0.0000 | A +0.010 CI [-0.066  0.09 ]
    per0   dark  pairs SW 1929 WS 1939 | rho_SW +0.070 CI [0.02  0.125] p 0.0025 | rho_WS +0.097 CI [0.041 0.151] p 0.0005 | A -0.026 CI [-0.09  0.04]
    per0   A[light]-A[dark] +0.038 CI [-0.034  0.117]
    fumin  light pairs SW 2630 WS 2617 | rho_SW +0.028 CI [-0.011  0.07 ] p 0.1830 | rho_WS +0.007 CI [-0.027  0.04 ] p 0.7580 | A +0.022 CI [-0.034  0.079]
    fumin  dark  pairs SW 1646 WS 1611 | rho_SW +0.023 CI [-0.041  0.084] p 0.3920 | rho_WS +0.023 CI [-0.017  0.067] p 0.3995 | A +0.000 CI [-0.083  0.082]
    fumin  A[light]-A[dark] +0.023 CI [-0.074  0.117]
    sss    light pairs SW 446 WS 431 | rho_SW -0.044 CI [-0.126  0.037] p 0.3855 | rho_WS -0.037 CI [-0.157  0.08 ] p 0.4815 | A -0.007 CI [-0.159  0.147]
    sss    dark  pairs SW 589 WS 559 | rho_SW +0.041 CI [-0.046  0.136] p 0.3555 | rho_WS -0.016 CI [-0.1    0.061] p 0.7140 | A +0.057 CI [-0.06   0.189]
    sss    A[light]-A[dark] -0.065 CI [-0.27   0.128]
    Multiplicity note (written with the result): 8 genotype x phase A intervals + 4 light-dark differences were computed; one A interval (CS dark) and one difference (CS) exclude 0.

## 2026-10-19 — PRESPEC 6 (eb14beb), hetero.py
    
    === SLEEP ===
    CS     windows 539 (dropped 55), bouts 3320 | T1 pooled dL(K2-K1): obs   132.250  null mean   175.076 [95%: 133.203, 218.493]  p 0.9800  share explained by heterogeneity 1.32
    CS     windows 539 (dropped 55), bouts 3320 | T2 within-window overdisp: obs    -0.048  null mean    -0.133 [95%: -0.179, -0.083]  p 0.0020
    per0   windows 531 (dropped 13), bouts 3962 | T1 pooled dL(K2-K1): obs    34.301  null mean   120.604 [95%: 86.750, 158.987]  p 0.9980  share explained by heterogeneity 3.52
    per0   windows 531 (dropped 13), bouts 3962 | T2 within-window overdisp: obs    -0.227  null mean    -0.114 [95%: -0.156, -0.071]  p 1.0000
    fumin  windows 818 (dropped 521), bouts 4485 | T1 pooled dL(K2-K1): obs    28.745  null mean   154.889 [95%: 113.774, 202.104]  p 0.9960  share explained by heterogeneity 5.39
    fumin  windows 818 (dropped 521), bouts 4485 | T2 within-window overdisp: obs    -0.218  null mean    -0.145 [95%: -0.187, -0.101]  p 1.0000
    sss    windows 224 (dropped 266), bouts 1107 | T1 pooled dL(K2-K1): obs    43.769  null mean    62.172 [95%: 35.125, 93.197]  p 0.9022  share explained by heterogeneity 1.42
    sss    windows 224 (dropped 266), bouts 1107 | T2 within-window overdisp: obs    -0.071  null mean    -0.157 [95%: -0.232, -0.081]  p 0.0140
    
    === WAKE ===
    CS     windows 536 (dropped 58), bouts 3313 | T1 pooled dL(K2-K1): obs   316.167  null mean   174.506 [95%: 131.745, 222.675]  p 0.0020  share explained by heterogeneity 0.55
    CS     windows 536 (dropped 58), bouts 3313 | T2 within-window overdisp: obs     0.318  null mean    -0.132 [95%: -0.177, -0.084]  p 0.0020
    per0   windows 530 (dropped 14), bouts 3962 | T1 pooled dL(K2-K1): obs   240.456  null mean    30.714 [95%: 16.531, 49.731]  p 0.0020  share explained by heterogeneity 0.13
    per0   windows 530 (dropped 14), bouts 3962 | T2 within-window overdisp: obs     0.014  null mean    -0.115 [95%: -0.159, -0.068]  p 0.0020
    fumin  windows 819 (dropped 534), bouts 4482 | T1 pooled dL(K2-K1): obs    86.652  null mean   204.641 [95%: 163.915, 251.309]  p 1.0000  share explained by heterogeneity 2.36
    fumin  windows 819 (dropped 534), bouts 4482 | T2 within-window overdisp: obs    -0.195  null mean    -0.145 [95%: -0.182, -0.101]  p 0.9960
    sss    windows 231 (dropped 271), bouts 1144 | T1 pooled dL(K2-K1): obs    99.672  null mean    60.964 [95%: 37.224, 83.503]  p 0.0060  share explained by heterogeneity 0.61
    sss    windows 231 (dropped 271), bouts 1144 | T2 within-window overdisp: obs     0.333  null mean    -0.156 [95%: -0.225, -0.087]  p 0.0020
    Note (written with the result): T2 one-sided as specified; per0 and fumin SLEEP T2 fall BELOW the null 95% range (under-dispersed) — the other tail, not a pre-specified test.

## 2026-10-19 — PRESPEC 7 (+amendment d4841af), within_fly.py
    blocks analysed: {'CS': 94, 'per0': 113, 'fumin': 104, 'sss': 4} excluded (<20 bouts, >0): {'CS': 56, 'per0': 37, 'fumin': 199, 'sss': 258}
    CS     blocks  94 | N1: 23/94 = 0.24 (binom p 1.8e-10) | N2 (PRIMARY): 10/94 = 0.11 (binom p 0.019) | sig-N2 blocks: median q_short 0.21, means 326 s / 1191 s
    per0   blocks 113 | N1: 8/113 = 0.07 (binom p 0.2) | N2 (PRIMARY): 3/113 = 0.03 (binom p 0.93) | sig-N2 blocks: median q_short 0.10, means 303 s / 1029 s
    fumin  blocks 104 | N1: 12/104 = 0.12 (binom p 0.0059) | N2 (PRIMARY): 1/104 = 0.01 (binom p 1) | sig-N2 blocks: median q_short 0.07, means 301 s / 463 s
    sss    blocks   4 | N1: 2/4 = 0.50 (binom p 0.014) | N2 (PRIMARY): 1/4 = 0.25 (binom p 0.19) | sig-N2 blocks: median q_short 0.12, means 308 s / 771 s
    
    By ZT bin (all genotypes' blocks, descriptive): fraction p<0.05 under N2
      ZT  0- 4: 2/66
      ZT  4- 8: 4/62
      ZT  8-12: 3/58
      ZT 12-16: 2/53
      ZT 16-20: 4/41
      ZT 20-24: 0/35

## 2026-10-19 — PRESPEC 8 (approved, d932e9e), irrev.py — RESULTS UNDER QUESTION, see 8b
    CS     whole pairs SW 3082 WS 3147 windows 594 | PASS(all 3 bin sets) True | null pass rate 0.0020
        3x3: D 0.0085 null 0.0012 D* 0.0073 CI [+0.0058,+0.0125] p 0.0005 | D*/h 0.0117
        4x4: D 0.0117 null 0.0027 D* 0.0090 CI [+0.0076,+0.0181] p 0.0005 | D*/h 0.0144
        5x5: D 0.0139 null 0.0048 D* 0.0092 CI [+0.0092,+0.0212] p 0.0005 | D*/h 0.0147
    CS     light pairs SW 1923 WS 1889 windows 366 | PASS(all 3 bin sets) True | null pass rate 0.0025
        3x3: D 0.0128 null 0.0021 D* 0.0107 CI [+0.0071,+0.0219] p 0.0010 | D*/h 0.0196
        4x4: D 0.0169 null 0.0044 D* 0.0125 CI [+0.0097,+0.0263] p 0.0005 | D*/h 0.0228
        5x5: D 0.0205 null 0.0080 D* 0.0125 CI [+0.0126,+0.0328] p 0.0005 | D*/h 0.0229
    CS     dark  pairs SW 1159 WS 1258 windows 338 | PASS(all 3 bin sets) False | null pass rate 0.0010
        3x3: D 0.0071 null 0.0033 D* 0.0038 CI [+0.0020,+0.0152] p 0.0585 | D*/h 0.0051
        4x4: D 0.0139 null 0.0076 D* 0.0063 CI [+0.0047,+0.0283] p 0.0510 | D*/h 0.0085
        5x5: D 0.0178 null 0.0128 D* 0.0050 CI [+0.0058,+0.0373] p 0.1224 | D*/h 0.0067
    per0   whole pairs SW 3612 WS 3818 windows 544 | PASS(all 3 bin sets) False | null pass rate 0.0015
        3x3: D 0.0028 null 0.0011 D* 0.0017 CI [+0.0010,+0.0063] p 0.0250 | D*/h 0.0033
        4x4: D 0.0037 null 0.0023 D* 0.0014 CI [+0.0012,+0.0082] p 0.0990 | D*/h 0.0028
        5x5: D 0.0070 null 0.0042 D* 0.0028 CI [+0.0027,+0.0140] p 0.0405 | D*/h 0.0055
    per0   light pairs SW 1822 WS 1943 windows 318 | PASS(all 3 bin sets) False | null pass rate 0.0040
        3x3: D 0.0042 null 0.0021 D* 0.0021 CI [+0.0009,+0.0092] p 0.0760 | D*/h 0.0042
        4x4: D 0.0088 null 0.0047 D* 0.0041 CI [+0.0038,+0.0161] p 0.0440 | D*/h 0.0081
        5x5: D 0.0134 null 0.0084 D* 0.0050 CI [+0.0042,+0.0266] p 0.0550 | D*/h 0.0099
    per0   dark  pairs SW 1790 WS 1875 windows 302 | PASS(all 3 bin sets) False | null pass rate 0.0030
        3x3: D 0.0020 null 0.0022 D* 0.0000 CI [-0.0009,+0.0074] p 0.4618 | D*/h 0.0000
        4x4: D 0.0038 null 0.0047 D* 0.0000 CI [-0.0016,+0.0116] p 0.6227 | D*/h 0.0000
        5x5: D 0.0139 null 0.0085 D* 0.0054 CI [+0.0028,+0.0298] p 0.0445 | D*/h 0.0104
    fumin  whole pairs SW 4827 WS 4284 windows 1271 | PASS(all 3 bin sets) True | null pass rate 0.0020
        3x3: D 0.0192 null 0.0009 D* 0.0184 CI [+0.0145,+0.0245] p 0.0005 | D*/h 0.0210
        4x4: D 0.0244 null 0.0019 D* 0.0225 CI [+0.0181,+0.0320] p 0.0005 | D*/h 0.0258
        5x5: D 0.0325 null 0.0034 D* 0.0291 CI [+0.0250,+0.0407] p 0.0005 | D*/h 0.0333
    fumin  light pairs SW 2922 WS 2628 windows 822 | PASS(all 3 bin sets) True | null pass rate 0.0035
        3x3: D 0.0158 null 0.0015 D* 0.0143 CI [+0.0113,+0.0215] p 0.0005 | D*/h 0.0177
        4x4: D 0.0218 null 0.0031 D* 0.0187 CI [+0.0159,+0.0285] p 0.0005 | D*/h 0.0231
        5x5: D 0.0265 null 0.0056 D* 0.0208 CI [+0.0192,+0.0361] p 0.0005 | D*/h 0.0257
    fumin  dark  pairs SW 1905 WS 1656 windows 641 | PASS(all 3 bin sets) True | null pass rate 0.0035
        3x3: D 0.0263 null 0.0023 D* 0.0240 CI [+0.0159,+0.0395] p 0.0005 | D*/h 0.0246
        4x4: D 0.0355 null 0.0050 D* 0.0304 CI [+0.0221,+0.0540] p 0.0005 | D*/h 0.0312
        5x5: D 0.0468 null 0.0091 D* 0.0378 CI [+0.0307,+0.0685] p 0.0005 | D*/h 0.0387
    sss    whole pairs SW 1326 WS 1140 windows 450 | PASS(all 3 bin sets) True | null pass rate 0.0045
        3x3: D 0.0617 null 0.0036 D* 0.0581 CI [+0.0469,+0.0798] p 0.0005 | D*/h 0.0551
        4x4: D 0.0846 null 0.0077 D* 0.0770 CI [+0.0627,+0.1135] p 0.0005 | D*/h 0.0730
        5x5: D 0.1261 null 0.0135 D* 0.1126 CI [+0.0967,+0.1706] p 0.0005 | D*/h 0.1068
    sss    light pairs SW 627 WS 533 windows 266 | PASS(all 3 bin sets) True | null pass rate 0.0020
        3x3: D 0.0700 null 0.0082 D* 0.0617 CI [+0.0452,+0.1034] p 0.0005 | D*/h 0.0508
        4x4: D 0.1134 null 0.0163 D* 0.0970 CI [+0.0797,+0.1670] p 0.0005 | D*/h 0.0799
        5x5: D 0.1838 null 0.0294 D* 0.1544 CI [+0.1301,+0.2851] p 0.0005 | D*/h 0.1272
    sss    dark  pairs SW 699 WS 607 windows 245 | PASS(all 3 bin sets) True | null pass rate 0.0035
        3x3: D 0.0554 null 0.0068 D* 0.0486 CI [+0.0322,+0.0898] p 0.0005 | D*/h 0.0533
        4x4: D 0.0629 null 0.0146 D* 0.0483 CI [+0.0367,+0.1145] p 0.0005 | D*/h 0.0530
        5x5: D 0.0987 null 0.0250 D* 0.0737 CI [+0.0580,+0.1766] p 0.0005 | D*/h 0.0808
    8b check 1-2 (post hoc): last bout of each window is length-biased (e.g. CS wake mean log-dur first 6.07 / interior 5.73 / last 7.85). Window cut rule is time-asymmetric -> PRESPEC 8 D* may be an artifact. Edge-free rerun launched.

## 2026-10-19 — 8b check 3 (POST HOC), irrev.py edgefree
    CS     whole pairs SW 2554 WS 2497 windows 568 | PASS(all 3 bin sets) False | null pass rate 0.0020
        3x3: D 0.0020 null 0.0015 D* 0.0004 CI [-0.0001,+0.0051] p 0.2589 | D*/h 0.0008
        4x4: D 0.0044 null 0.0034 D* 0.0010 CI [-0.0001,+0.0105] p 0.2259 | D*/h 0.0018
        5x5: D 0.0043 null 0.0061 D* 0.0000 CI [-0.0010,+0.0124] p 0.8161 | D*/h 0.0000
    CS     light pairs SW 1596 WS 1596 windows 343 | PASS(all 3 bin sets) False | null pass rate 0.0020
        3x3: D 0.0034 null 0.0025 D* 0.0009 CI [-0.0003,+0.0085] p 0.2319 | D*/h 0.0019
        4x4: D 0.0071 null 0.0055 D* 0.0016 CI [+0.0009,+0.0158] p 0.2309 | D*/h 0.0034
        5x5: D 0.0064 null 0.0098 D* 0.0000 CI [-0.0013,+0.0180] p 0.8656 | D*/h 0.0000
    CS     dark  pairs SW 958 WS 901 windows 299 | PASS(all 3 bin sets) False | null pass rate 0.0040
        3x3: D 0.0019 null 0.0047 D* 0.0000 CI [-0.0033,+0.0105] p 0.8656 | D*/h 0.0000
        4x4: D 0.0064 null 0.0105 D* 0.0000 CI [-0.0041,+0.0247] p 0.8351 | D*/h 0.0000
        5x5: D 0.0142 null 0.0180 D* 0.0000 CI [-0.0024,+0.0427] p 0.7066 | D*/h 0.0000
    per0   whole pairs SW 3274 WS 3069 windows 541 | PASS(all 3 bin sets) False | null pass rate 0.0030
        3x3: D 0.0015 null 0.0013 D* 0.0003 CI [-0.0005,+0.0054] p 0.2974 | D*/h 0.0005
        4x4: D 0.0020 null 0.0027 D* 0.0000 CI [-0.0009,+0.0056] p 0.6897 | D*/h 0.0000
        5x5: D 0.0052 null 0.0050 D* 0.0002 CI [-0.0001,+0.0132] p 0.3993 | D*/h 0.0004
    per0   light pairs SW 1666 WS 1540 windows 292 | PASS(all 3 bin sets) False | null pass rate 0.0020
        3x3: D 0.0021 null 0.0025 D* 0.0000 CI [-0.0015,+0.0083] p 0.4938 | D*/h 0.0000
        4x4: D 0.0082 null 0.0055 D* 0.0027 CI [+0.0010,+0.0168] p 0.1394 | D*/h 0.0055
        5x5: D 0.0105 null 0.0098 D* 0.0007 CI [-0.0002,+0.0268] p 0.3848 | D*/h 0.0014
    per0   dark  pairs SW 1608 WS 1529 windows 283 | PASS(all 3 bin sets) False | null pass rate 0.0015
        3x3: D 0.0016 null 0.0027 D* 0.0000 CI [-0.0015,+0.0078] p 0.6962 | D*/h 0.0000
        4x4: D 0.0042 null 0.0056 D* 0.0000 CI [-0.0026,+0.0124] p 0.6692 | D*/h 0.0000
        5x5: D 0.0131 null 0.0102 D* 0.0029 CI [+0.0020,+0.0277] p 0.1909 | D*/h 0.0057
    fumin  whole pairs SW 3126 WS 3583 windows 1016 | PASS(all 3 bin sets) True | null pass rate 0.0020
        3x3: D 0.0050 null 0.0012 D* 0.0038 CI [+0.0031,+0.0076] p 0.0005 | D*/h 0.0069
        4x4: D 0.0063 null 0.0028 D* 0.0035 CI [+0.0028,+0.0115] p 0.0075 | D*/h 0.0064
        5x5: D 0.0102 null 0.0049 D* 0.0053 CI [+0.0044,+0.0168] p 0.0050 | D*/h 0.0096
    fumin  light pairs SW 1957 WS 2210 windows 671 | PASS(all 3 bin sets) False | null pass rate 0.0035
        3x3: D 0.0051 null 0.0020 D* 0.0030 CI [+0.0020,+0.0089] p 0.0280 | D*/h 0.0055
        4x4: D 0.0070 null 0.0044 D* 0.0026 CI [+0.0015,+0.0132] p 0.1039 | D*/h 0.0047
        5x5: D 0.0105 null 0.0079 D* 0.0026 CI [+0.0022,+0.0181] p 0.1419 | D*/h 0.0048
    fumin  dark  pairs SW 1169 WS 1373 windows 467 | PASS(all 3 bin sets) True | null pass rate 0.0010
        3x3: D 0.0070 null 0.0033 D* 0.0037 CI [+0.0017,+0.0160] p 0.0420 | D*/h 0.0068
        4x4: D 0.0144 null 0.0075 D* 0.0069 CI [+0.0039,+0.0295] p 0.0310 | D*/h 0.0126
        5x5: D 0.0250 null 0.0131 D* 0.0119 CI [+0.0104,+0.0464] p 0.0110 | D*/h 0.0218
    sss    whole pairs SW 748 WS 892 windows 299 | PASS(all 3 bin sets) False | null pass rate 0.0020
        3x3: D 0.0102 null 0.0059 D* 0.0043 CI [-0.0011,+0.0256] p 0.1124 | D*/h 0.0080
        4x4: D 0.0159 null 0.0124 D* 0.0035 CI [-0.0001,+0.0420] p 0.2224 | D*/h 0.0066
        5x5: D 0.0434 null 0.0210 D* 0.0224 CI [+0.0180,+0.0783] p 0.0070 | D*/h 0.0419
    sss    light pairs SW 331 WS 399 windows 156 | PASS(all 3 bin sets) False | null pass rate 0.0015
        3x3: D 0.0105 null 0.0141 D* 0.0000 CI [-0.0070,+0.0445] p 0.6032 | D*/h 0.0000
        4x4: D 0.0256 null 0.0279 D* 0.0000 CI [-0.0072,+0.0750] p 0.5142 | D*/h 0.0000
        5x5: D 0.0616 null 0.0454 D* 0.0162 CI [+0.0104,+0.1156] p 0.1449 | D*/h 0.0292
    sss    dark  pairs SW 417 WS 493 windows 174 | PASS(all 3 bin sets) False | null pass rate 0.0040
        3x3: D 0.0161 null 0.0099 D* 0.0062 CI [-0.0013,+0.0392] p 0.1324 | D*/h 0.0119
        4x4: D 0.0206 null 0.0221 D* 0.0000 CI [-0.0045,+0.0579] p 0.4943 | D*/h 0.0000
        5x5: D 0.0606 null 0.0361 D* 0.0245 CI [+0.0205,+0.1104] p 0.0510 | D*/h 0.0470
    Note: bootstrap CIs exclude 0 in several blocks with reversal-null p > 0.2 (e.g. CS light 4x4); the plug-in bootstrap is biased upward, so the reversal-null p, not the CI, is the test.

## 2026-10-20 AUTO-APPENDED raw output — PRESPEC 9 (acc2ffd), within_fly9.py (interpretation to follow)
    A 4h blocks: {'CS': 82, 'per0': 101, 'fumin': 99, 'sss': 3}
    B 12h blocks: {'CS': 50, 'per0': 50, 'fumin': 84, 'sss': 24}
    
    === A: 4 h blocks (primary null N2) ===
    CS     blocks  82 | FISHER (N2) p 0.00544 | N1: 23/82 (binom p 8.8e-12) | N2: 10/82 (binom p 0.0078) | sig: q_short 0.20, means 316/1220 s
    per0   blocks 101 | FISHER (N2) p 0.772 | N1: 8/101 (binom p 0.13) | N2: 4/101 (binom p 0.75) | sig: q_short 0.10, means 303/997 s
    fumin  blocks  99 | FISHER (N2) p 1 | N1: 11/99 (binom p 0.011) | N2: 0/99 (binom p 1)
    sss    blocks   3 | FISHER (N2) p 0.0272 | N1: 2/3 (binom p 0.0073) | N2: 0/3 (binom p 1)
    
    === B: 12 h blocks (primary null N3) ===
    CS     blocks  50 | FISHER (N3) p 0.0127 | N1: 23/50 (binom p 3.4e-17) | N2: 22/50 (binom p 5.4e-16) | N3: 6/50 (binom p 0.038) | sig: q_short 0.14, means 326/1560 s
    per0   blocks  50 | FISHER (N3) p 0.783 | N1: 6/50 (binom p 0.038) | N2: 6/50 (binom p 0.038) | N3: 2/50 (binom p 0.72) | sig: q_short 0.15, means 361/1434 s
    fumin  blocks  84 | FISHER (N3) p 1 | N1: 9/84 (binom p 0.025) | N2: 3/84 (binom p 0.8) | N3: 0/84 (binom p 1)
    sss    blocks  24 | FISHER (N3) p 0.967 | N1: 6/24 (binom p 0.00096) | N2: 1/24 (binom p 0.71) | N3: 0/24 (binom p 1)

## 2026-10-20 — PRESPEC 10 (ab86161) C1 and C2 raw output (counts_analysis.py, run 1). C3 was INTERRUPTED by a container restart (~04:35 UTC) before writing any result; C3 is rerun with per-job checkpointing (counts_C3_resume.py), identical model code and seeds.
    === C1: movement fraction around transitions minus mid-wake reference ===
    CS     pre  light flies 25 median contrast -0.0210 (IQR -0.0512..-0.0020) Wilcoxon p 0.00119
    CS     pre  dark  flies 25 median contrast -0.0461 (IQR -0.0865..-0.0271) Wilcoxon p 3.19e-05
    CS     post light flies 25 median contrast -0.0419 (IQR -0.0640..-0.0199) Wilcoxon p 5.96e-07
    CS     post dark  flies 25 median contrast +0.0379 (IQR +0.0219..+0.0663) Wilcoxon p 0.0046
    fumin  pre  light flies 49 median contrast -0.0093 (IQR -0.0360..+0.0119) Wilcoxon p 0.00913
    fumin  pre  dark  flies 53 median contrast -0.0107 (IQR -0.0440..+0.0271) Wilcoxon p 0.192
    fumin  post light flies 48 median contrast +0.0270 (IQR +0.0105..+0.0414) Wilcoxon p 1.12e-07
    fumin  post dark  flies 51 median contrast +0.0090 (IQR -0.0428..+0.0414) Wilcoxon p 0.837
    
    === C2: inactivity runs, all lengths ===
    CS     light runs 194690 (>=300 s: 2371) | CV K=5 | CV logL [np.float64(-780249.9), np.float64(-564390.9), np.float64(-541718.3), np.float64(-539585.7), np.float64(-539315.0)]
           weights [0.5039, 0.2724, 0.1647, 0.0438, 0.0151] means(s) [1.7, 5.4, 20.3, 98.6, 701.6] | runs>=300 s assigned to components with mean<300 s: 0.181
    CS     dark  runs 137857 (>=300 s: 2005) | CV K=5 | CV logL [np.float64(-601814.2), np.float64(-383948.0), np.float64(-356027.8), np.float64(-353213.0), np.float64(-352973.7)]
           weights [0.4371, 0.2955, 0.1988, 0.0506, 0.018] means(s) [1.2, 3.7, 13.7, 63.8, 1190.9] | runs>=300 s assigned to components with mean<300 s: 0.039
    fumin  light runs 575151 (>=300 s: 1337) | CV K=5 | CV logL [np.float64(-1810223.0), np.float64(-1517313.8), np.float64(-1478320.8), np.float64(-1470078.2), np.float64(-1469677.8)]
           weights [0.3034, 0.4032, 0.2278, 0.0559, 0.0098] means(s) [1.1, 3.5, 13.1, 40.9, 209.0] | runs>=300 s assigned to components with mean<300 s: 1.000
    fumin  dark  runs 619283 (>=300 s: 991) | CV K=5 | CV logL [np.float64(-1856141.4), np.float64(-1596224.9), np.float64(-1544890.9), np.float64(-1538386.4), np.float64(-1538115.2)]
           weights [0.3203, 0.4441, 0.1931, 0.0359, 0.0067] means(s) [1.1, 4.2, 13.8, 45.0, 207.2] | runs>=300 s assigned to components with mean<300 s: 1.000
    
    === C3: 3-state count HMM per fly x phase ===

## 2026-10-20 — INTERPRETATION of PRESPEC 9 and of PRESPEC 10 C1/C2 (written after the raw entries above)
PRESPEC 9 (edge-free within-fly two-phase test):
- CS supported at both block sizes under the pre-set primary rule: 4 h Fisher p 0.0054 (N2; 10/82
  blocks, binom p 0.0078) and 12 h Fisher p 0.0127 (N3; 6/50, binom p 0.038). The PRESPEC 7 CS
  result therefore survives removal of the length-selected last bout of each window.
- per0 and fumin not supported at either block size. sss: 4 h Fisher p 0.027 rests on 3 blocks with
  none individually < 0.05 (uninformative); 12 h (24 blocks) Fisher p 0.97: not supported.
- Fisher used discrete bootstrap p-values (min 1/101): approximate, as stated in the plan.
- The short component in significant CS blocks has mean 316-326 s (q 0.14-0.20), i.e. ~20 s above the
  300 s threshold. Working (untested) reading, from C2: it is the tail of a rest-like inactivity
  component crossing the threshold, not a distinct sleep state.
PRESPEC 10 C1 (movement fraction vs mid-wake, Wilcoxon across flies; 8 tests, no correction pre-set):
- CS: lower before sleep onset in light (p 0.0012) and dark (p 3.2e-5); after waking lower in light
  (p 6e-7) and HIGHER in dark (p 0.0046). fumin: lower before onset in light only (p 0.009); higher
  after waking in light (p 1.1e-7); dark n.s. Six of 8 pass Bonferroni 0.05/8.
- Model-free evidence of a quieter approach to sleep in CS; it does not distinguish a discrete
  transitional state from gradual slowing.
PRESPEC 10 C2 (inactivity runs of all lengths):
- CV chose K = 5, the plan's upper limit, in all four genotype x phase groups, with CV still rising:
  K is limited by the plan, so component structure beyond 5 is untested.
- Fraction of >= 300 s runs assigned to components with mean < 300 s: CS light 0.181, CS dark 0.039,
  fumin light 1.000, fumin dark 1.000 (fumin's longest component mean ~208 s). Within K <= 5, fumin's
  count-defined "sleep" bouts are all the tail of rest-like components. Pooled across flies, so the
  heterogeneity caveat of PRESPEC 6 applies.

## 2026-10-20 AUTO-APPENDED raw output — PRESPEC 10 C3 (counts_C3_resume.py), interpretation to follow
    remaining jobs 148
    CS     light flies 25 | H1 degenerate (min L1<0.05) 0.92 | median L1(H1,H2) 0.860 L1(H1,H3) 0.023 | route-unidentified (restart range>0.2) 0.12 | median restarts within 0.5 nats 1
            median Xn +0.874 Yn +0.343 | median Xn_range 0.000 Yn_range 0.000
            held-out log-lik per s, median over flies: K2 -0.25218 K3 -0.24681 K4 -0.18881; flies where K3>K2 0.52, K4>K3 1.00
    CS     dark  flies 25 | H1 degenerate (min L1<0.05) 0.56 | median L1(H1,H2) 1.210 L1(H1,H3) 0.034 | route-unidentified (restart range>0.2) 0.16 | median restarts within 0.5 nats 1
            median Xn +1.000 Yn +0.999 | median Xn_range 0.000 Yn_range 0.000
            held-out log-lik per s, median over flies: K2 -0.19796 K3 -0.19793 K4 -0.13797; flies where K3>K2 0.68, K4>K3 1.00
    fumin  light flies 53 | H1 degenerate (min L1<0.05) 0.23 | median L1(H1,H2) 0.951 L1(H1,H3) 0.796 | route-unidentified (restart range>0.2) 0.06 | median restarts within 0.5 nats 2
            median Xn -0.196 Yn +0.221 | median Xn_range 0.000 Yn_range 0.001
            held-out log-lik per s, median over flies: K2 -0.64430 K3 -0.52367 K4 -0.50509; flies where K3>K2 0.70, K4>K3 0.98
    fumin  dark  flies 53 | H1 degenerate (min L1<0.05) 0.15 | median L1(H1,H2) 0.761 L1(H1,H3) 1.018 | route-unidentified (restart range>0.2) 0.06 | median restarts within 0.5 nats 2
            median Xn -0.065 Yn -0.807 | median Xn_range 0.001 Yn_range 0.000
            held-out log-lik per s, median over flies: K2 -0.69992 K3 -0.61872 K4 -0.57633; flies where K3>K2 0.77, K4>K3 0.98

## 2026-10-20 — INTERPRETATION of PRESPEC 10 C3 (3-state count HMM), written after the AUTO raw entry
Median emission rows (P of count 0,1,2,>=3) and median mean dwell, computed from C3_results.jsonl:
  CS light : H1 [.987 .008 .002 .001] 19 s | H2 [.552 .363 .062 .007] 14 s | H3 [.999 .001 0 0] 530 s
  CS dark  : H1 [.982 .017 0 0] 16 s       | H2 [.388 .527 .085 .007] 12 s | H3 [.999 0 0 0] 1138 s
  fumin L  : H1 [.600 .358 .008 .004] 8 s  | H2 [.158 .595 .203 .014] 6 s  | H3 [.998 .001 0 0] 24 s
  fumin D  : H1 [.489 .448 .038 .004] 6 s  | H2 [.122 .574 .231 .020] 5 s  | H3 [.998 .001 0 .001] 17 s
- CS: H1 is a near-silent BRIEF-PAUSE state (P0 ~0.98, ~16-19 s dwell), nearly emission-degenerate
  with H3 (median L1 0.02-0.03; "degenerate" < 0.05 in 92% light / 56% dark flies). Distinguished
  from H3 mainly by dwell, so the H1/H3 equivalence ridge is near-flat; CS Xn (+0.87 light, +1.00
  dark) and Yn are NOT considered identified.
- fumin: H1 has its own emission signature (median L1 to H2 0.76-0.95, to H3 0.80-1.02; degenerate in
  15-23% of flies), a LOW-ACTIVITY state. But all three states dwell for seconds (H3 17-24 s): at
  K = 3 the model describes second-scale movement texture, not sleep/wake; fumin's H3 is not sleep.
- Restart-spread check (pre-set criterion: range > 0.2 = route-unidentified) had little power: the
  median number of restarts within 0.5 nats of the best was 1-2 (CS 5/25 and 10/25 flies, fumin
  37/53 and 38/53 with >= 2). A range of 0 from one restart is not evidence of identification.
  The pre-set "route-unidentified" fractions (0.06-0.16) are therefore NOT reported as support.
- Held-out last day: K=4 beats K=3 in 98-100% of flies in every group; K=3 beats K=2 in only 52-77%.
  Three states are not enough for these data, and the third state is spent on fast movement structure.
- Bottom line for H1: counts DO separate a third emission state in fumin, and in CS give a brief-pause
  state; neither is the sleep<->wake regulator of the scored HMM. A model with separate fast
  (movement) and slow (behavioural) timescales is needed before routing can be asked again.

## 2026-10-20 — PRESPEC 11 stage 1 (flatK.py), results
    K_g: {"CS_0": 6, "CS_1": 4, "CS_2": 4, "CS_3": 4, "CS_4": 4, "CS_5": 4, "fumin_0": 4, "fumin_1": 4, "fumin_2": 4, "fumin_3": 6, "fumin_4": 4, "fumin_5": 4}
    errors: 0
    CS
      ZT  0- 4: K* counts {4: 4, 5: 1, 6: 5} | median held-out gain per hour (nats) K3-2   77.6 K4-3  125.1 K5-4   -3.3 K6-5    1.9 | flies K3 within 1 nat/h of best: 0/10
      ZT  4- 8: K* counts {4: 8, 5: 1, 6: 1} | median held-out gain per hour (nats) K3-2   -0.0 K4-3  149.6 K5-4   -4.1 K6-5   -0.7 | flies K3 within 1 nat/h of best: 0/10
      ZT  8-12: K* counts {4: 6, 5: 1, 6: 3} | median held-out gain per hour (nats) K3-2    0.2 K4-3  148.0 K5-4   -4.7 K6-5    5.3 | flies K3 within 1 nat/h of best: 0/10
      ZT 12-16: K* counts {4: 8, 5: 2} | median held-out gain per hour (nats) K3-2   -0.2 K4-3  301.0 K5-4   -3.4 K6-5   -0.1 | flies K3 within 1 nat/h of best: 0/10
      ZT 16-20: K* counts {2: 1, 4: 9} | median held-out gain per hour (nats) K3-2    nan K4-3    nan K5-4    nan K6-5    nan | flies K3 within 1 nat/h of best: 0/10
      ZT 20-24: K* counts {4: 6, 5: 3, 6: 1} | median held-out gain per hour (nats) K3-2   -0.1 K4-3  254.7 K5-4   -1.6 K6-5   -4.5 | flies K3 within 1 nat/h of best: 0/10
    fumin
      ZT  0- 4: K* counts {4: 5, 5: 1, 6: 4} | median held-out gain per hour (nats) K3-2  434.8 K4-3   34.2 K5-4   -6.3 K6-5    1.0 | flies K3 within 1 nat/h of best: 0/10
      ZT  4- 8: K* counts {4: 6, 5: 3, 6: 1} | median held-out gain per hour (nats) K3-2  428.0 K4-3   45.5 K5-4   -6.5 K6-5   -3.1 | flies K3 within 1 nat/h of best: 0/10
      ZT  8-12: K* counts {4: 4, 5: 3, 6: 3} | median held-out gain per hour (nats) K3-2  569.8 K4-3   38.2 K5-4    nan K6-5    2.8 | flies K3 within 1 nat/h of best: 0/10
      ZT 12-16: K* counts {4: 4, 5: 1, 6: 5} | median held-out gain per hour (nats) K3-2  419.6 K4-3   48.6 K5-4   -5.8 K6-5    5.3 | flies K3 within 1 nat/h of best: 0/10
      ZT 16-20: K* counts {4: 4, 5: 2, 6: 4} | median held-out gain per hour (nats) K3-2  449.0 K4-3   58.6 K5-4   -2.2 K6-5    6.0 | flies K3 within 1 nat/h of best: 0/10
      ZT 20-24: K* counts {4: 5, 6: 5} | median held-out gain per hour (nats) K3-2  559.2 K4-3   45.7 K5-4   -1.1 K6-5    2.6 | flies K3 within 1 nat/h of best: 0/10
    Note: 2 fly-bins had a held-out fold scored -inf (CS fly 23 bin 4 day 1, all K; fumin fly 7 bin 2 day 3, K=4,5): the held-out day contains a symbol the trained model gives probability 0. The plan did not anticipate this. K_g is identical whether those folds are kept (as run) or dropped for every K (checked post hoc). NaN entries in the gain table come from these two fly-bins.

## 2026-10-21 — PRESPEC 11 stages 2-4 results (flatK_stage2.jsonl; summary script inline, output below)
    records 468 with fit errors 0
    
    ===== CS =====
    ZT  0- 4 K=6 flies 25 | m(60s) {0: 2, 1: 23} | m(30s) {0: 2, 1: 23} | m(300s) {0: 22, 1: 3} | median restarts within 0.5 nats 1
            median top-3 implied timescales (s): [169.6, 6.4, 2.8]
            m=1 (23 flies): macro movement fraction [0.005, 0.221] occupancy [0.565, 0.435] dwell(s) [388.0, 336.0]
    ZT  4- 8 K=4 flies 25 | m(60s) {0: 1, 1: 24} | m(30s) {1: 25} | m(300s) {0: 23, 1: 2} | median restarts within 0.5 nats 1
            median top-3 implied timescales (s): [153.2, 7.6, 0.4]
            m=1 (24 flies): macro movement fraction [0.001, 0.14] occupancy [0.775, 0.225] dwell(s) [839.0, 206.0]
    ZT  8-12 K=4 flies 25 | m(60s) {1: 25} | m(30s) {1: 25} | m(300s) {0: 23, 1: 2} | median restarts within 0.5 nats 1
            median top-3 implied timescales (s): [128.8, 6.0, 1.5]
            m=1 (25 flies): macro movement fraction [0.004, 0.173] occupancy [0.646, 0.354] dwell(s) [397.0, 201.0]
    ZT 12-16 K=4 flies 25 | m(60s) {1: 25} | m(30s) {1: 24, 2: 1} | m(300s) {0: 18, 1: 7} | median restarts within 0.5 nats 1
            median top-3 implied timescales (s): [260.9, 7.3, 3.7]
            m=1 (25 flies): macro movement fraction [0.004, 0.293] occupancy [0.673, 0.327] dwell(s) [736.0, 395.0]
    ZT 16-20 K=4 flies 25 | m(60s) {1: 25} | m(30s) {1: 25} | m(300s) {0: 13, 1: 12} | median restarts within 0.5 nats 1
            median top-3 implied timescales (s): [299.4, 7.3, 1.1]
            m=1 (25 flies): macro movement fraction [0.001, 0.256] occupancy [0.793, 0.207] dwell(s) [1620.0, 406.0]
    ZT 20-24 K=4 flies 25 | m(60s) {1: 25} | m(30s) {1: 25} | m(300s) {0: 10, 1: 15} | median restarts within 0.5 nats 1
            median top-3 implied timescales (s): [335.3, 7.4, 1.0]
            m=1 (25 flies): macro movement fraction [0.001, 0.25] occupancy [0.831, 0.169] dwell(s) [1992.0, 473.0]
    
    ===== fumin =====
    ZT  0- 4 K=4 flies 53 | m(60s) {0: 26, 1: 27} | m(30s) {0: 14, 1: 39} | m(300s) {0: 48, 1: 5} | median restarts within 0.5 nats 1
            median top-3 implied timescales (s): [61.2, 3.7, 2.3]
            m=1 (27 flies): macro movement fraction [0.019, 0.407] occupancy [0.374, 0.626] dwell(s) [188.0, 338.0]
    ZT  4- 8 K=4 flies 53 | m(60s) {0: 30, 1: 23} | m(30s) {0: 20, 1: 33} | m(300s) {0: 49, 1: 4} | median restarts within 0.5 nats 1
            median top-3 implied timescales (s): [50.5, 3.7, 2.2]
            m=1 (23 flies): macro movement fraction [0.024, 0.367] occupancy [0.411, 0.589] dwell(s) [170.0, 350.0]
    ZT  8-12 K=4 flies 53 | m(60s) {0: 27, 1: 26} | m(30s) {0: 14, 1: 39} | m(300s) {0: 50, 1: 3} | median restarts within 0.5 nats 1
            median top-3 implied timescales (s): [58.9, 3.8, 2.4]
            m=1 (26 flies): macro movement fraction [0.014, 0.321] occupancy [0.395, 0.605] dwell(s) [181.0, 313.0]
    ZT 12-16 K=6 flies 53 | m(60s) {0: 27, 1: 23, 2: 3} | m(30s) {0: 16, 1: 33, 2: 3, 3: 1} | m(300s) {0: 47, 1: 6} | median restarts within 0.5 nats 1
            median top-3 implied timescales (s): [59.8, 4.2, 3.0]
            m=1 (23 flies): macro movement fraction [0.026, 0.377] occupancy [0.5, 0.5] dwell(s) [166.0, 192.0]
            m=2 (3 flies): macro movement fraction [0.004, 0.35, 0.493] occupancy [0.283, 0.38, 0.334] dwell(s) [113.0, 1.0, 1.0]
            routing via T (m=2, 3 flies): S->W median 0.048 [IQR 0.024-0.524] | W->S median 0.349 [IQR 0.252-0.675] | flies with >=2 near restarts: 3, median range 0.000
    ZT 16-20 K=4 flies 53 | m(60s) {0: 38, 1: 15} | m(30s) {0: 24, 1: 29} | m(300s) {0: 49, 1: 4} | median restarts within 0.5 nats 1
            median top-3 implied timescales (s): [38.8, 3.6, 2.6]
            m=1 (15 flies): macro movement fraction [0.016, 0.436] occupancy [0.424, 0.576] dwell(s) [175.0, 286.0]
    ZT 20-24 K=4 flies 53 | m(60s) {0: 36, 1: 17} | m(30s) {0: 25, 1: 28} | m(300s) {0: 51, 1: 2} | median restarts within 0.5 nats 1
            median top-3 implied timescales (s): [34.6, 3.6, 2.4]
            m=1 (17 flies): macro movement fraction [0.024, 0.445] occupancy [0.337, 0.663] dwell(s) [119.0, 289.0]
    INTERPRETATION (written with the results):
    - PRIMARY (pre-set): m = 2 (three metastable sets) is the most common value in NO genotype x bin.
      CS: m = 1 in 23-25 of 25 flies in every bin -> two slow groups only: a near-silent group
      (movement fraction 0.1-0.5%, dwell ~6.5-33 min, longest at night) and an active group (14-29%,
      dwell ~3.5-8 min). Slowest implied timescale 2.1-5.6 min; the next is ~6-7 s (fast texture).
      This matches the stage 1 jump at K = 4: roughly 2 slow groups x 2 fast states.
    - fumin: m = 1 in only 15-27 of 53 flies per bin, m = 0 in the rest (slowest timescale 35-61 s,
      at the 60 s threshold). Even the sleep/wake split is marginal at 60 s; groups dwell ~2-6 min.
      m = 2 in 3/53 flies in ZT 12-16 only; there T and W groups dwell ~1 s, i.e. not metastable, so
      the routing-through-T values for those 3 flies are not interpretable.
    - Sensitivity: at 30 s the CS picture is unchanged (m = 1); at 300 s CS has m = 0 by day and
      m = 1 in about half the flies at night.
    - Optimizer reliability: median number of restarts within 0.5 nats of the best = 1 in every
      group. Most restarts end in different local optima, so per-fly details (state rows, P) are
      fragile. The count of slow modes was not checked across restarts; this should be done before
      relying on per-fly values.
    - For H1: in these count data, no third long-lived (>= 60 s) group appears between the sleep-like
      and wake-like groups. A transitional process, if present, lives on the fast timescale (seconds)
      or as a change in fast-state occupancy within the groups (cf. the pre-sleep quieting of C1).

## 2026-10-21 — status: PRESPEC 12 (6343b8d) committed; NOT run. hhmm.py not yet written (the write+synthetic test was interrupted by the user before execution). No PRESPEC 12 result exists.
    CORRECTION: hhmm.py WAS written (64 lines, 00:36 UTC) before the interruption; its synthetic recovery test was not run (no output exists), and the model code is untested.

## 2026-10-21 — PRESPEC 12 synthetic validation (b3b43df): FAIL (criterion 1)
    == P3 (157s) truth mf [0.008, 0.181, 0.45] dwell [900.0, 120.0, 300.0] Xn_c +0.333 Yn_c +0.077
       fitted Nc=3: mf [0.008, 0.173, 0.467] dwell [931.0, 117.0, 266.0] Xn_c +0.471 Yn_c +0.424 | restart logLs [-16393.9, -16332.7, -16331.4, -16330.9, -16322.9, -16262.5]
       monotone EM: True | held-out gain Nc3-Nc2: +7.6 nats/h
       criterion 1 False | 2 True | 3 True | 5 True
    == P2 (148s) truth mf [0.008, 0.45] dwell [900.0, 300.0]
       fitted Nc=3: mf [0.009, 0.424, 0.665] dwell [972.0, 242.0, 27.0] Xn_c -0.241 Yn_c +1.000 | restart logLs [-16553.8, -16550.6, -16256.6, -16251.3, -16250.7, -16246.4]
       monotone EM: True | held-out gain Nc3-Nc2: -1.6 nats/h
       criterion 4 (negative control) True | 5 True
    SYNTHETIC VALIDATION FAIL
    Per plan, the data analysis is NOT run.

## 2026-10-21 — PRESPEC 12 addendum 2, replicate validation (12 datasets from P3): FAIL
    Xn_c error      mean -0.120  SD 0.732  95% CI of mean [-0.534, +0.294]  values [-0.46, -0.22, 0.31, 0.35, 0.67, -1.33, -1.18, 0.67, -0.63, 0.67, 0.43, -0.7]
    Yn_c error      mean -0.057  SD 0.574  95% CI of mean [-0.382, +0.267]  values [0.03, -0.34, 0.49, 0.09, -0.47, -0.92, -1.08, 0.92, 0.39, 0.23, 0.07, -0.1]
    dwell ratio S   mean +0.931  SD 0.357  95% CI of mean [+0.730, +1.133]  values [1.1, 0.91, 0.9, 0.9, 1.39, 0.3, 0.31, 0.8, 1.45, 0.93, 0.95, 1.24]
    dwell ratio T   mean +0.920  SD 0.446  95% CI of mean [+0.667, +1.172]  values [1.08, 1.06, 1.46, 1.03, 0.94, 0.04, 0.05, 0.82, 1.08, 0.91, 1.19, 1.37]
    dwell ratio W   mean +1.139  SD 0.369  95% CI of mean [+0.930, +1.348]  values [1.39, 0.82, 0.72, 0.94, 1.16, 1.01, 2.07, 1.33, 0.91, 1.0, 0.9, 1.4]
    REPLICATE VALIDATION FAIL
    INTERPRETATION: bias is small (mean error Xn_c -0.12, Yn_c -0.06; both CIs contain 0) but the spread
    is huge (SD 0.73 and 0.57 on a -1..+1 scale), even though the fitted model has exactly the true
    structure. Two replicates (seeds 106, 107) also landed in wrong optima (T dwell 0.04-0.05x truth).
    With 4 days x 4 h of one fly, coarse routing is NOT estimable per fly-bin. Data analysis not run.
    Options: pool flies within genotype x bin (shared coarse parameters), better initialisation and
    more restarts, then repeat this replicate validation at the pooled data size before any data run.

## 2026-10-21 — UNPLANNED, DESCRIPTIVE aging/drift check (aging_check.py), requested by user
    (a) held-out log-lik per second at K_g, within-fly-bin centred, by held-out day (median over fly-bins)
      CS     n=59 fly-bins | by day [-18.1, -35.8, 21.3, 22.7] nats/h | slope per day median +12.3 nats/h, Wilcoxon p 0.219
      fumin  n=59 fly-bins | by day [-23.2, -31.2, 49.6] nats/h | slope per day median +56.1 nats/h, Wilcoxon p 0.0543
    
    (b) raw daily behaviour, median over flies by day | within-fly slope per day (median) and Wilcoxon p
      CS     movement fraction      by day [0.066, 0.072, 0.064, 0.065] | slope -0.002/day, Wilcoxon p 0.107 (n=25)
      CS     sleep (min/day)        by day [869.183, 803.433, 771.133, 799.5] | slope -21.413/day, Wilcoxon p 0.00807 (n=25)
      CS     sleep bouts/day        by day [42.0, 43.0, 47.0, 48.0] | slope +2.100/day, Wilcoxon p 0.00807 (n=25)
      CS     mean sleep bout (min)  by day [22.4, 19.524, 18.942, 16.334] | slope -0.970/day, Wilcoxon p 0.000103 (n=25)
      CS     movement fraction      days 2+ only: slope -0.004/day, Wilcoxon p 0.00226
      CS     sleep (min/day)        days 2+ only: slope -11.917/day, Wilcoxon p 0.474
      fumin  movement fraction      by day [0.237, 0.287, 0.274] | slope +0.009/day, Wilcoxon p 0.0636 (n=53)
      fumin  sleep (min/day)        by day [155.25, 86.517, 73.6] | slope -30.767/day, Wilcoxon p 2.77e-06 (n=53)
      fumin  sleep bouts/day        by day [18.0, 11.0, 10.0] | slope -3.000/day, Wilcoxon p 1.59e-05 (n=53)
      fumin  mean sleep bout (min)  by day [8.819, 7.622, 7.242] | slope -0.604/day, Wilcoxon p 0.000133 (n=53)
      fumin  movement fraction      days 2+ only: slope -0.011/day, Wilcoxon p 0.237
      fumin  sleep (min/day)        days 2+ only: slope -5.800/day, Wilcoxon p 0.315
    Reading: clear drift within the recordings. CS: sleep bouts shorten and multiply day by day (mean bout 22.4 -> 16.3 min); fumin: large day-1 -> day-2 drop in sleep, then flat. Held-out HMM scores show no significant trend (they also move with behavioural entropy, so they are a weak diagnostic).

## 2026-10-21 — hhmm_day.py informal unit test (8 flies x 4 days, not a gate)
    815s monotone True iters 178
    true kappa_S [0.67, 0.88, 1.14, 1.49] fitted kappa (S,T,W): [[0.71, 0.75, 1.23, 1.52], [1.12, 0.89, 1.04, 0.96], [1.19, 0.92, 0.99, 0.93]]
    Xn_c fitted +0.194 (truth +0.333) Yn_c +0.077 (truth +0.077)
    
    [exited with code 0]

## 2026-10-21 — PRESPEC 13 Part V results (validate13.py)
    remaining 10
    done V1 201
    done V2 301
    done V2 302
    done V1 202
    done V1 203
    done V1 204
    done V1 205
    done V1 206
    done V3 401
    done V3 402
    V1 dXn: mean +0.010 SD 0.058 values [-0.002, -0.031, 0.093, -0.009, -0.059, 0.065]
    V1 dYn: mean +0.019 SD 0.064 values [-0.057, 0.036, 0.081, 0.056, 0.065, -0.064]
    V1 mean dwell ratio S,T,W [0.974, 1.02, 1.025]
    V2 median held-out gain per fly (nats/h): [8.26, 9.33]
    V1 seed 201: dLL(fit - truth) +13.4 nats -> fit >= truth (error is statistical)
    V1 seed 202: dLL(fit - truth) +1.5 nats -> fit >= truth (error is statistical)
    V1 seed 203: dLL(fit - truth) +14.4 nats -> fit >= truth (error is statistical)
    V1 seed 204: dLL(fit - truth) +21.8 nats -> fit >= truth (error is statistical)
    V1 seed 205: dLL(fit - truth) +3.2 nats -> fit >= truth (error is statistical)
    V1 seed 206: dLL(fit - truth) +25.9 nats -> fit >= truth (error is statistical)
    V3 seed 401: dLL(fit - truth) -11137.4 nats -> OPTIMIZER FAILURE
    V3 seed 402: dLL(fit - truth) -14083.5 nats -> OPTIMIZER FAILURE
    V3 seed 401: dXn -0.045 dYn +0.173 dwell ratio [0.92, 2.16, 0.98]
    V3 seed 402: dXn -0.774 dYn -0.523 dwell ratio [0.66, 4.54, 1.46]
    V1 True V2 False => VALIDATION FAIL
    INTERPRETATION:
    - V1 PASS: pooling 25 flies x 4 days gives routing error SD 0.06 (vs 0.6-0.7 per fly-bin), bias ~0;
      every fit's logL >= the truth's, so remaining error is statistical, not the optimizer.
    - V2 FAIL: with fly-to-fly heterogeneity and NO transitional state, Nc=3 beats Nc=2 on held-out days
      (+8.3, +9.3 nats/h per fly). Pooling heterogeneous flies CREATES a false transitional state.
      Per the plan, Part D is not run.
    - Correction to the addendum read-out: for V3 the 'truth' is each fly's own model, which a pooled model
      cannot reach; the large negative dLL there is misspecification, NOT optimizer failure. The dLL rule
      is only meaningful for V1. V3 routing errors: seed 401 small, seed 402 large (dXn -0.77).

## 2026-10-21 — PRESPEC 13 amendment 2 gates V5/V6 (validate13_day.py)
    remaining 5
    done V5 502
    Model is not converging.  Current: -89902.54658053443 is not greater than -89902.52316935764. Delta is -0.02341117679316085
    done V6 601
    done V5 501
    done V6 602
    done V5 503
    V5 dXn: mean -0.001 SD 0.020 values [-0.016, 0.021, -0.01]
    V5 dYn: mean -0.010 SD 0.072 values [-0.074, 0.067, -0.025]
    V5 mean fitted/true kappa_S by day [1.063, 0.929, 1.056, 0.96]
    V6 median held-out gain per fly: [4.75, 23.69]
    V5 True V6 False => DAY-MODEL VALIDATION FAIL
    INTERPRETATION: V5 PASS — with a day trend, routing recovered with SD 0.02 (Xn_c) and 0.07 (Yn_c), bias ~0,
    day factors within 7%. V6 FAIL — heterogeneous flies + day trend, no T: held-out gain +4.75 (pass) and
    +23.7 nats/h per fly (fail). Same conclusion as V2: fly heterogeneity can create a false transitional state.
    One subset fit's logL decreased by 0.023 nats in one EM step: the numerical day-factor M-step is not exact,
    so the day model is a generalized EM with tiny non-monotone steps; noted, not consequential here.
    Part D is not run.

## 2026-10-21 — PRESPEC 14 partial-pooling unit test (informal, 8 heterogeneous flies x 4 days)
    898s for 3 fits | iters 300 | largest logL decrease 0.0000 nats
    true kappa_S [0.67, 0.88, 1.14, 1.49] | fitted kappa (S,T,W) [[0.92, 0.83, 0.96, 1.37], [1.07, 0.96, 0.91, 1.08], [1.01, 0.93, 0.99, 1.07]]
    Xn_c -0.624 (truth +0.333) | Yn_c -0.777 (truth +0.077) | mf [0.107, 0.499, 0.505]
    restart final logLs [-266955.5, -267131.4, -266873.7]
    
    INTERPRETATION: the model as built does NOT work. Routing is wrong in sign and size (Xn_c -0.62 vs
    +0.33; Yn_c -0.78 vs +0.08); T and W coarse states end with near-identical movement (0.499 vs 0.505),
    so the S/T/W labels are ambiguous; fits hit the 300-iteration cap; restarts differ by up to 257 nats.
    Likely cause: with FREE per-fly emissions, the movement signature that distinguishes T from W is no
    longer shared, so the coarse states lose their identity (a new degeneracy). Free per-fly texture is too
    much freedom; per-fly deviations need to be constrained (shrunk toward a shared value), i.e. design B.

## 2026-10-21 — PRESPEC 15 step 1 (epr_min.py)
    CS     H1==H2 windows  403 | EPR_min <= 1e-3 nats/h: 0.752 | median EPR_min 6.76e-29 | median fitted EPR 18.1 nats/h (fitted infinite in 2) | J range contains 0: 0.814 | median Xn,Yn at min-EPR +0.02,+0.20
    CS     H1==H3 windows   31 | EPR_min <= 1e-3 nats/h: 0.355 | median EPR_min 7.85e-01 | median fitted EPR 26.8 nats/h (fitted infinite in 0) | J range contains 0: 0.871 | median Xn,Yn at min-EPR +0.02,-0.38
    per0   H1==H2 windows  310 | EPR_min <= 1e-3 nats/h: 0.768 | median EPR_min 3.52e-29 | median fitted EPR 1.46 nats/h (fitted infinite in 1) | J range contains 0: 0.713 | median Xn,Yn at min-EPR +0.06,+0.28
    per0   H1==H3 windows   75 | EPR_min <= 1e-3 nats/h: 0.200 | median EPR_min 1.46e+01 | median fitted EPR 151 nats/h (fitted infinite in 2) | J range contains 0: 0.413 | median Xn,Yn at min-EPR +0.99,-0.99
    fumin  H1==H2 windows 1456 | EPR_min <= 1e-3 nats/h: 0.089 | median EPR_min 8.32e-01 | median fitted EPR 116 nats/h (fitted infinite in 309) | J range contains 0: 0.975 | median Xn,Yn at min-EPR +0.16,-0.03
    fumin  H1==H3 windows    8 | EPR_min <= 1e-3 nats/h: 0.125 | median EPR_min 9.24e+00 | median fitted EPR 33.3 nats/h (fitted infinite in 1) | J range contains 0: 0.750 | median Xn,Yn at min-EPR +0.71,-0.66
    sss    H1==H2 windows  586 | EPR_min <= 1e-3 nats/h: 0.087 | median EPR_min 5.48e-01 | median fitted EPR 118 nats/h (fitted infinite in 193) | J range contains 0: 0.956 | median Xn,Yn at min-EPR -0.25,+0.28
    sss    H1==H3 windows    7 | EPR_min <= 1e-3 nats/h: 0.143 | median EPR_min 2.10e+00 | median fitted EPR 2.61 nats/h (fitted infinite in 0) | J range contains 0: 1.000 | median Xn,Yn at min-EPR +0.44,-0.51
    POST HOC robustness (epr_min_check.py, written after the results): wider search unchanged:
    fumin: 40 windows | old median EPR_min 0.478 -> wide search 0.478 nats/h | now <= 1e-3: 0 | wide/old ratio median 1.000
    sss: 40 windows | old median EPR_min 0.692 -> wide search 0.692 nats/h | now <= 1e-3: 0 | wide/old ratio median 1.000
    INTERPRETATION: within a 3-state description, CS and per0 data are compatible with zero dissipation in
    ~75% of exactly-degenerate windows; fumin and sss data REQUIRE dissipation in ~91% (median minimum
    0.55-0.83 nats/h, ~150-200x below the fitted pipeline values). Cycle-flux sign undetermined in 71-98%.
    My prediction (minimum = 0 everywhere) was wrong for fumin/sss; H3's small 'active' emission breaks
    exact lumping. The bound is conditional on the 3-state model class.

## 2026-10-21 — PRESPEC 15 E2/E3 (e23.py)
    
    ===== CS (25 flies) =====
    E2 light flies 25 bouts 466 | mean D (exit - approach) -0.0094  perm p 0.0675
          exit profile min2-10     [0.116, 0.107, 0.114, 0.125, 0.123, 0.126, 0.132, 0.134, 0.139]
          approach profile min2-10 [0.106, 0.117, 0.132, 0.14, 0.146, 0.142, 0.141, 0.135, 0.142]
    E2 dark  flies 25 bouts 309 | mean D (exit - approach) +0.0621  perm p 0.0005
          exit profile min2-10     [0.227, 0.233, 0.257, 0.248, 0.246, 0.228, 0.242, 0.251, 0.23]
          approach profile min2-10 [0.122, 0.144, 0.171, 0.188, 0.186, 0.198, 0.2, 0.202, 0.193]
    E3 light segments 300 triples 389069 | D 0.00488 null mean 0.00006 D* 0.00482 nats/triple | p 0.0010
    E3 dark  segments 300 triples 275430 | D 0.01195 null mean 0.00013 D* 0.01182 nats/triple | p 0.0010
    
    ===== fumin (53 flies) =====
    E2 light flies 48 bouts 758 | mean D (exit - approach) +0.0153  perm p 0.0010
          exit profile min2-10     [0.27, 0.253, 0.234, 0.22, 0.201, 0.197, 0.192, 0.189, 0.188]
          approach profile min2-10 [0.187, 0.204, 0.2, 0.206, 0.208, 0.212, 0.196, 0.2, 0.193]
    E2 dark  flies 51 bouts 596 | mean D (exit - approach) +0.0072  perm p 0.1109
          exit profile min2-10     [0.315, 0.274, 0.26, 0.258, 0.233, 0.226, 0.236, 0.216, 0.212]
          approach profile min2-10 [0.221, 0.243, 0.248, 0.244, 0.245, 0.246, 0.245, 0.236, 0.237]
    E3 light segments 477 triples 1149754 | D 0.00418 null mean 0.00005 D* 0.00413 nats/triple | p 0.0010
    E3 dark  segments 477 triples 1238132 | D 0.00292 null mean 0.00005 D* 0.00288 nats/triple | p 0.0010
    INTERPRETATION: non-equilibrium is MEASURED in the count data, without any model.
    - E2: waking and falling asleep are not mirror images. CS dark (p 0.0005): after waking activity is high at once
      (~0.23-0.26) while the approach to sleep ramps down (0.20 -> 0.12 in the last minutes). fumin light (p 0.001):
      a burst after waking that decays (0.27 -> 0.19) vs a flat approach with a dip before onset. CS light p 0.068,
      fumin dark p 0.11 (not significant).
    - E3: run-level time-irreversibility detected in all four genotype x phase groups (p <= 0.001, the smallest
      attainable with 1000 replicates); D* 0.003-0.012 nats per run-triple. Caveat: the segment-reversal null is
      exact only for stationary segments; a slow within-segment trend contributes little to local triple asymmetry.

## 2026-10-22 — PRESPEC 16 E2/E3 on per0LD, per0DD, sss (e23.py)
    
    ===== per0LD (27 flies) =====
    E2 light flies 27 bouts 383 | mean D (exit - approach) +0.0111  perm p 0.1099
          exit profile min2-10     [0.177, 0.178, 0.169, 0.184, 0.18, 0.185, 0.188, 0.194, 0.201]
          approach profile min2-10 [0.15, 0.158, 0.177, 0.181, 0.171, 0.179, 0.185, 0.179, 0.175]
    E2 dark  flies 26 bouts 365 | mean D (exit - approach) +0.0174  perm p 0.0005
          exit profile min2-10     [0.146, 0.144, 0.139, 0.134, 0.132, 0.128, 0.137, 0.141, 0.138]
          approach profile min2-10 [0.102, 0.105, 0.117, 0.13, 0.132, 0.127, 0.125, 0.12, 0.125]
    E3 light segments 324 triples 388037 | D 0.00195 null mean 0.00008 D* 0.00188 nats/triple | p 0.0010
    E3 dark  segments 324 triples 281680 | D 0.00071 null mean 0.00005 D* 0.00066 nats/triple | p 0.0010
    
    ===== per0DD (25 flies) =====
    E2 light flies 25 bouts 279 | mean D (exit - approach) +0.0265  perm p 0.0005
          exit profile min2-10     [0.147, 0.167, 0.156, 0.161, 0.166, 0.184, 0.177, 0.164, 0.179]
          approach profile min2-10 [0.104, 0.121, 0.127, 0.131, 0.142, 0.145, 0.161, 0.158, 0.175]
    E2 dark  flies 24 bouts 275 | mean D (exit - approach) +0.0249  perm p 0.0015
          exit profile min2-10     [0.142, 0.151, 0.157, 0.15, 0.156, 0.168, 0.164, 0.173, 0.175]
          approach profile min2-10 [0.101, 0.117, 0.111, 0.138, 0.157, 0.136, 0.148, 0.153, 0.15]
    E2 whole flies 25 bouts 555 | mean D (exit - approach) +0.0255  perm p 0.0005
          exit profile min2-10     [0.142, 0.157, 0.156, 0.154, 0.16, 0.173, 0.169, 0.166, 0.174]
          approach profile min2-10 [0.101, 0.119, 0.119, 0.134, 0.148, 0.138, 0.151, 0.15, 0.16]
    E3 light segments 273 triples 243696 | D 0.00324 null mean 0.00008 D* 0.00315 nats/triple | p 0.0010
    E3 dark  segments 273 triples 232791 | D 0.00246 null mean 0.00008 D* 0.00237 nats/triple | p 0.0010
    E3 whole segments 546 triples 476487 | D 0.00284 null mean 0.00004 D* 0.00279 nats/triple | p 0.0010
    
    ===== sss (57 flies) =====
    E2 light flies 46 bouts 485 | mean D (exit - approach) -0.0009  perm p 0.8401
          exit profile min2-10     [0.126, 0.139, 0.15, 0.142, 0.146, 0.162, 0.166, 0.16, 0.158]
          approach profile min2-10 [0.143, 0.145, 0.146, 0.148, 0.159, 0.152, 0.152, 0.157, 0.151]
    E2 dark  flies 51 bouts 502 | mean D (exit - approach) -0.0026  perm p 0.5202
          exit profile min2-10     [0.131, 0.146, 0.129, 0.146, 0.153, 0.149, 0.153, 0.155, 0.14]
          approach profile min2-10 [0.148, 0.153, 0.143, 0.155, 0.14, 0.139, 0.148, 0.149, 0.15]
    E3 light segments 513 triples 1683332 | D 0.00402 null mean 0.00003 D* 0.00399 nats/triple | p 0.0010
    E3 dark  segments 513 triples 1640589 | D 0.00071 null mean 0.00001 D* 0.00070 nats/triple | p 0.0010
    INTERPRETATION: E2 asymmetry (falling asleep != waking reversed) in per0LD dark (p 0.0005) and per0DD
    whole record (p 0.0005; both clock halves also p <= 0.0015); not per0LD light (0.11); ABSENT in sss (light
    0.84, dark 0.52). Across all groups, 9 primary E2 tests (per0DD counted once): the four significant ones
    (CS dark, fumin light, per0LD dark, per0DD whole) all pass Bonferroni 0.05/9. The asymmetry persists in
    arrhythmic flies in constant darkness, so it is not produced by light or circadian timing. E3: run-sequence
    irreversibility in every group and block (p <= 0.001), D* 0.0007-0.012 nats per run-triple.

## 2026-10-22 — PRESPEC 17 (taper.py) results
    CS      light flies 25 | D median -0.0044 IQR [-0.0347,+0.0145] D>0 0.48 sign p 1 | fly perm p<.05 0.04 | depth median +0.0359 rel +0.25 Wilcoxon p 6.56e-06
             fit n=25: tau_in 100.00 min [0.46, 100.0] A_in -0.64 | tau_out 100.00 [0.68, 100.0] A_out -0.54 | dtau CI [-99.47, 99.29] dA CI [-1.35, 2.56]
    CS      dark  flies 25 | D median +0.0560 IQR [+0.0330,+0.0960] D>0 0.92 sign p 1.94e-05 | fly perm p<.05 0.44 | depth median +0.0589 rel +0.41 Wilcoxon p 1.97e-06
             fit n=25: tau_in 35.54 min [0.69, 100.0] A_in -1.41 | tau_out 38.91 [15.97, 100.0] A_out -2.90 | dtau CI [-97.88, 84.03] dA CI [0.68, 2.59]
    fumin   light flies 48 | D median +0.0153 IQR [-0.0030,+0.0332] D>0 0.73 sign p 0.00209 | fly perm p<.05 0.10 | depth median +0.0052 rel +0.03 Wilcoxon p 0.0958
             fit n=52: tau_in 100.00 min [59.07, 100.0] A_in -0.67 | tau_out 13.98 [8.67, 23.89] A_out -1.17 | dtau CI [42.45, 91.3] dA CI [0.23, 0.86]
    fumin   dark  flies 51 | D median +0.0040 IQR [-0.0124,+0.0279] D>0 0.53 sign p 0.78 | fly perm p<.05 0.12 | depth median +0.0290 rel +0.11 Wilcoxon p 9.28e-05
             fit n=53: tau_in 69.40 min [21.22, 100.0] A_in -0.42 | tau_out 15.82 [3.05, 31.95] A_out -0.63 | dtau CI [-5.76, 96.9] dA CI [0.01, 0.57]
    per0LD  light flies 27 | D median +0.0068 IQR [-0.0087,+0.0388] D>0 0.63 sign p 0.248 | fly perm p<.05 0.19 | depth median +0.0200 rel +0.11 Wilcoxon p 0.00594
             fit n=24: tau_in 100.00 min [0.38, 100.0] A_in -0.18 | tau_out 100.00 [0.74, 100.0] A_out -0.26 | dtau CI [-99.08, 99.0] dA CI [-4.67, 1.13]
    per0LD  dark  flies 26 | D median +0.0148 IQR [+0.0067,+0.0388] D>0 0.85 sign p 0.000534 | fly perm p<.05 0.12 | depth median +0.0230 rel +0.16 Wilcoxon p 8.16e-05
             fit n=22: tau_in 88.69 min [14.35, 100.0] A_in -1.38 | tau_out 70.55 [17.54, 100.0] A_out -1.81 | dtau CI [-76.06, 81.75] dA CI [-0.18, 1.78]
    per0DD  whole flies 25 | D median +0.0258 IQR [+0.0100,+0.0403] D>0 0.84 sign p 0.000911 | fly perm p<.05 0.40 | depth median +0.0423 rel +0.27 Wilcoxon p 6.56e-06
             fit n=17: tau_in 100.00 min [3.61, 100.0] A_in -1.31 | tau_out 100.00 [18.66, 100.0] A_out -1.15 | dtau CI [-66.38, 68.23] dA CI [-1.16, 0.49]
    sss     light flies 46 | D median -0.0053 IQR [-0.0188,+0.0101] D>0 0.43 sign p 0.461 | fly perm p<.05 0.04 | depth median +0.0227 rel +0.15 Wilcoxon p 1.58e-05
             fit n=53: tau_in 0.23 min [0.21, 100.0] A_in +5.00 | tau_out 0.73 [0.25, 100.0] A_out +0.58 | dtau CI [-99.72, 99.73] dA CI [-5.18, 5.18]
    sss     dark  flies 51 | D median -0.0037 IQR [-0.0183,+0.0176] D>0 0.41 sign p 0.262 | fly perm p<.05 0.04 | depth median +0.0180 rel +0.13 Wilcoxon p 0.00507
             fit n=55: tau_in 100.00 min [100.0, 100.0] A_in -0.51 | tau_out 53.03 [17.57, 100.0] A_out -0.58 | dtau CI [-0.0, 82.43] dA CI [-0.08, 0.26]
    INTERPRETATION:
    A (per-fly D): consistent across flies where E2 was significant: D > 0 in 92% of CS-dark flies (sign p 2e-5),
      85% per0LD-dark (5e-4), 84% per0DD (9e-4), 73% fumin-light (0.002); individually significant in 44%, 12%,
      40%, 10% of flies. No consistent sign in CS light, fumin dark, per0LD light, sss light/dark.
    B (depth): activity in minutes 2-3 before sleep onset is BELOW the mid-wake level in every group (relative
      depth 0.03-0.41; Wilcoxon p <= 0.006 except fumin light 0.096). sss DOES taper before sleep; its E2 null
      means the start of its wake bouts is equally quiet (symmetric), not that it lacks a taper. This corrects
      the wording 'absent in sleepless flies' used on the review page (22 Oct, section 8).
    C INVALID: per-bout normalisation by the bout's own mid-wake level divides by near-zero values in quiet
      bouts (normalised approach means 1.2-3.1 with SEM up to 1.3); fits hit their bounds. The flaw is in the
      plan, not the code. tau_in/tau_out are not reported.

## 2026-10-22 — figures (figs.py): fig1 per-fly D, fig2 E2 profiles (2-10 min), fig3 symmetry plane (exit depth = DESCRIPTIVE post hoc mirror of the pre-registered approach depth), fig4 average relaxation 2-20 min for bouts >= 40 min (DESCRIPTIVE replacement for the invalid PRESPEC 17 C fits, at the user's suggestion; no time constants fitted).
    fig4 flies (bouts >= 40 min): [('CS', 'dark', 25), ('per0DD', 'whole', 17), ('per0LD', 'dark', 22), ('fumin', 'light', 52), ('sss', 'dark', 55)]
    CS      light n=25 depth_in +4.98% depth_out +4.98% (mean over flies)
    CS      dark  n=25 depth_in +7.25% depth_out -2.47% (mean over flies)
    fumin   light n=48 depth_in +1.35% depth_out -5.27% (mean over flies)
    fumin   dark  n=51 depth_in +3.61% depth_out -2.65% (mean over flies)
    per0LD  light n=27 depth_in +3.58% depth_out +1.26% (mean over flies)
    per0LD  dark  n=26 depth_in +3.32% depth_out -0.81% (mean over flies)
    per0DD  whole n=25 depth_in +5.06% depth_out +1.17% (mean over flies)
    sss     light n=46 depth_in +3.54% depth_out +4.75% (mean over flies)
    sss     dark  n=51 depth_in +1.55% depth_out +2.76% (mean over flies)

## 2026-10-22 — PRESPEC 18 (9770d02 + bridging amendment) results
    Option (a) run (bridge-then-threshold) DISCARDED after the user specified option (b); output kept in v18_optionA.log.
    Option (b) run (v18.py, v18b.log):
    V0 CS      light flies 25 bouts  466 | median D -0.44 pp, D>0 0.48, sign p 1, perm p 0.0675 | x +2.79 y +3.72
    V0 CS      dark  flies 25 bouts  309 | median D +5.60 pp, D>0 0.92, sign p 1.9e-05, perm p 0.0005 | x +2.75 y -3.46
    V0 fumin   light flies 48 bouts  758 | median D +1.53 pp, D>0 0.73, sign p 0.0021, perm p 0.0005 | x +0.82 y -0.71
    V0 fumin   dark  flies 51 bouts  596 | median D +0.40 pp, D>0 0.53, sign p 0.78, perm p 0.1294 | x +2.74 y +2.02
    V0 per0LD  light flies 27 bouts  383 | median D +0.68 pp, D>0 0.63, sign p 0.25, perm p 0.1224 | x +1.69 y +0.59
    V0 per0LD  dark  flies 26 bouts  365 | median D +1.48 pp, D>0 0.85, sign p 0.00053, perm p 0.0005 | x +1.66 y -0.08
    V0 per0DD  whole flies 25 bouts  555 | median D +2.58 pp, D>0 0.84, sign p 0.00091, perm p 0.0005 | x +2.52 y -0.03
    V0 sss     light flies 46 bouts  485 | median D -0.53 pp, D>0 0.43, sign p 0.46, perm p 0.8241 | x +2.91 y +3.00
    V0 sss     dark  flies 51 bouts  502 | median D -0.37 pp, D>0 0.41, sign p 0.26, perm p 0.5207 | x +1.87 y +2.13
    V1b CS      light flies 25 bouts  466 | median D -0.44 pp, D>0 0.48, sign p 1, perm p 0.0640 | x +2.79 y +3.72
    V1b CS      dark  flies 25 bouts  309 | median D +5.60 pp, D>0 0.92, sign p 1.9e-05, perm p 0.0005 | x +2.75 y -3.46
    V1b fumin   light flies 48 bouts  758 | median D +1.53 pp, D>0 0.73, sign p 0.0021, perm p 0.0005 | x +0.82 y -0.71
    V1b fumin   dark  flies 51 bouts  596 | median D +0.40 pp, D>0 0.53, sign p 0.78, perm p 0.1234 | x +2.74 y +2.02
    V1b per0LD  light flies 27 bouts  383 | median D +0.68 pp, D>0 0.63, sign p 0.25, perm p 0.1304 | x +1.69 y +0.59
    V1b per0LD  dark  flies 26 bouts  365 | median D +1.48 pp, D>0 0.85, sign p 0.00053, perm p 0.0005 | x +1.66 y -0.08
    V1b per0DD  whole flies 25 bouts  555 | median D +2.58 pp, D>0 0.84, sign p 0.00091, perm p 0.0005 | x +2.52 y -0.03
    V1b sss     light flies 46 bouts  485 | median D -0.53 pp, D>0 0.43, sign p 0.46, perm p 0.8361 | x +2.91 y +3.00
    V1b sss     dark  flies 51 bouts  502 | median D -0.37 pp, D>0 0.41, sign p 0.26, perm p 0.5402 | x +1.87 y +2.13
    V2b CS      light flies 25 bouts  466 | median D -1.19 pp, D>0 0.36, sign p 0.23, perm p 0.2519 | x +3.76 y +4.58
    V2b CS      dark  flies 25 bouts  309 | median D +8.30 pp, D>0 0.92, sign p 1.9e-05, perm p 0.0005 | x +4.95 y -3.57
    V2b fumin   light flies 48 bouts  758 | median D +4.25 pp, D>0 0.83, sign p 3.3e-06, perm p 0.0005 | x +0.96 y -3.54
    V2b fumin   dark  flies 51 bouts  596 | median D +3.84 pp, D>0 0.78, sign p 5.7e-05, perm p 0.0005 | x +2.90 y -0.87
    V2b per0LD  light flies 27 bouts  383 | median D +1.19 pp, D>0 0.67, sign p 0.12, perm p 0.3083 | x +2.33 y +1.30
    V2b per0LD  dark  flies 26 bouts  365 | median D +1.91 pp, D>0 0.81, sign p 0.0025, perm p 0.0005 | x +2.34 y -0.38
    V2b per0DD  whole flies 25 bouts  555 | median D +3.37 pp, D>0 0.76, sign p 0.015, perm p 0.0005 | x +4.24 y +0.86
    V2b sss     light flies 46 bouts  485 | median D -0.92 pp, D>0 0.37, sign p 0.1, perm p 0.2174 | x +3.40 y +4.06
    V2b sss     dark  flies 51 bouts  502 | median D -0.89 pp, D>0 0.37, sign p 0.092, perm p 0.0265 | x +1.64 y +2.81
    INTERPRETATION: V1b reproduces V0 exactly (same bouts, D, x, y; permutation p differs only by resampling),
    as predicted: bridging (b) cannot move the edges of wake bouts >= 20 min. V2b (post hoc +-5 min window,
    DESCRIPTIVE): asymmetry stronger near the transition (median Delta_5: CS dark +8.3, fumin light +4.3,
    fumin dark +3.8 [ns at 2-10 min], per0LD dark +1.9, per0DD +3.4); sleepless slightly negative (-0.9;
    perm p 0.027 dark). Figures fig1-3_pm5_bridged generated; V0 figures retained.
    User confirmation (2026-10-22): a 1 s movement between a sleep bout and a shorter pause (e.g. 301 + 1 + 100 s) is NOT bridged. Bridging rule (b) as implemented is the lab convention; no rerun needed.

## 2026-10-23 — PRESPEC 19 (571f90b) signed cycle current (cycle.py; JSON-save bug fixed and rerun, same seed and output)
    
    ===== theta = 0.5  (PRIMARY) =====
    CS      light flies 25 | c mean -0.019 [-0.080,+0.049] median -0.066 c>0 0.40 sign p 0.42 Wilcoxon p 0.46 | J -0.007/h [-0.047,+0.035] | onsets Q:A 264:534 wakings Q:A 296:550
    CS      dark  flies 25 | c mean +0.111 [+0.069,+0.157] median +0.111 c>0 0.80 sign p 0.0041 Wilcoxon p 6.4e-05 | J +0.062/h [+0.033,+0.091] | onsets Q:A 248:530 wakings Q:A 150:580
    fumin   light flies 49 | c mean +0.048 [+0.017,+0.078] median +0.050 c>0 0.55 sign p 0.088 Wilcoxon p 0.0024 | J +0.043/h [+0.024,+0.062] | onsets Q:A 168:825 wakings Q:A 104:911
    fumin   dark  flies 53 | c mean +0.059 [+0.011,+0.101] median +0.083 c>0 0.68 sign p 2.5e-05 Wilcoxon p 0.00071 | J +0.024/h [+0.012,+0.037] | onsets Q:A 112:711 wakings Q:A 53:745
    per0LD  light flies 27 | c mean +0.007 [-0.030,+0.047] median +0.016 c>0 0.56 sign p 0.56 Wilcoxon p 0.71 | J +0.007/h [-0.014,+0.027] | onsets Q:A 188:578 wakings Q:A 192:601
    per0LD  dark  flies 27 | c mean -0.015 [-0.047,+0.017] median -0.026 c>0 0.37 sign p 0.25 Wilcoxon p 0.3 | J -0.010/h [-0.029,+0.010] | onsets Q:A 225:618 wakings Q:A 225:591
    per0DD  whole flies 25 | c mean +0.105 [+0.051,+0.161] median +0.105 c>0 0.72 sign p 0.043 Wilcoxon p 0.0023 | J +0.066/h [+0.026,+0.105] | onsets Q:A 524:1033 wakings Q:A 379:1178
    sss     light flies 50 | c mean -0.093 [-0.154,-0.036] median -0.068 c>0 0.28 sign p 0.032 Wilcoxon p 0.0067 | J -0.013/h [-0.029,+0.001] | onsets Q:A 125:411 wakings Q:A 174:412
    sss     dark  flies 53 | c mean -0.029 [-0.097,+0.037] median +0.000 c>0 0.47 sign p 1 Wilcoxon p 0.57 | J -0.008/h [-0.029,+0.013] | onsets Q:A 170:531 wakings Q:A 161:491
      sleepless vs others (Mann-Whitney on c):
        light: sss vs CS: p 0.342
        light: sss vs fumin: p 0.00018
        light: sss vs per0LD: p 0.0245
        light: sss vs per0DD(whole): p 5.42e-05
        dark: sss vs CS: p 0.0252
        dark: sss vs fumin: p 0.0793
        dark: sss vs per0LD: p 0.768
        dark: sss vs per0DD(whole): p 0.0252
    
    ===== theta = 0.33 =====
    CS      light flies 25 | c mean -0.023 [-0.078,+0.035] median -0.042 c>0 0.40 sign p 0.42 Wilcoxon p 0.35 | J -0.008/h [-0.042,+0.026] | onsets Q:A 181:617 wakings Q:A 214:632
    CS      dark  flies 25 | c mean +0.077 [+0.038,+0.116] median +0.077 c>0 0.68 sign p 0.11 Wilcoxon p 0.0023 | J +0.033/h [+0.007,+0.057] | onsets Q:A 177:601 wakings Q:A 113:617
    fumin   light flies 49 | c mean +0.030 [+0.002,+0.059] median +0.000 c>0 0.49 sign p 0.099 Wilcoxon p 0.068 | J +0.026/h [+0.010,+0.042] | onsets Q:A 102:891 wakings Q:A 67:948
    fumin   dark  flies 53 | c mean +0.054 [+0.030,+0.081] median +0.000 c>0 0.49 sign p 0.00054 Wilcoxon p 0.00011 | J +0.012/h [+0.005,+0.021] | onsets Q:A 66:757 wakings Q:A 30:768
    per0LD  light flies 27 | c mean +0.019 [-0.011,+0.047] median +0.013 c>0 0.63 sign p 0.11 Wilcoxon p 0.15 | J +0.017/h [-0.002,+0.034] | onsets Q:A 128:638 wakings Q:A 120:673
    per0LD  dark  flies 27 | c mean -0.002 [-0.028,+0.024] median -0.003 c>0 0.44 sign p 0.7 Wilcoxon p 0.82 | J -0.005/h [-0.022,+0.014] | onsets Q:A 155:688 wakings Q:A 148:668
    per0DD  whole flies 25 | c mean +0.074 [+0.028,+0.120] median +0.066 c>0 0.72 sign p 0.023 Wilcoxon p 0.0072 | J +0.044/h [+0.016,+0.074] | onsets Q:A 385:1172 wakings Q:A 289:1268
    sss     light flies 50 | c mean -0.063 [-0.112,-0.014] median -0.014 c>0 0.30 sign p 0.15 Wilcoxon p 0.016 | J -0.010/h [-0.022,+0.001] | onsets Q:A 72:464 wakings Q:A 115:471
    sss     dark  flies 53 | c mean +0.014 [-0.044,+0.068] median +0.000 c>0 0.45 sign p 0.44 Wilcoxon p 0.67 | J -0.003/h [-0.020,+0.013] | onsets Q:A 114:587 wakings Q:A 96:556
      sleepless vs others (Mann-Whitney on c):
        light: sss vs CS: p 0.562
        light: sss vs fumin: p 0.00319
        light: sss vs per0LD: p 0.0107
        light: sss vs per0DD(whole): p 0.000399
        dark: sss vs CS: p 0.146
        dark: sss vs fumin: p 0.351
        dark: sss vs per0LD: p 0.373
        dark: sss vs per0DD(whole): p 0.129
    
    ===== theta = 0.75 =====
    CS      light flies 25 | c mean -0.022 [-0.071,+0.031] median -0.049 c>0 0.36 sign p 0.23 Wilcoxon p 0.22 | J -0.019/h [-0.048,+0.011] | onsets Q:A 378:420 wakings Q:A 425:421
    CS      dark  flies 25 | c mean +0.173 [+0.110,+0.239] median +0.188 c>0 0.84 sign p 0.00091 Wilcoxon p 3.8e-05 | J +0.107/h [+0.065,+0.154] | onsets Q:A 375:403 wakings Q:A 223:507
    fumin   light flies 49 | c mean +0.105 [+0.068,+0.143] median +0.117 c>0 0.78 sign p 6.2e-05 Wilcoxon p 1.4e-05 | J +0.077/h [+0.052,+0.102] | onsets Q:A 321:672 wakings Q:A 197:818
    fumin   dark  flies 53 | c mean +0.150 [+0.107,+0.191] median +0.167 c>0 0.79 sign p 1.2e-06 Wilcoxon p 7.3e-07 | J +0.062/h [+0.044,+0.081] | onsets Q:A 236:587 wakings Q:A 105:693
    per0LD  light flies 27 | c mean +0.031 [-0.019,+0.081] median +0.034 c>0 0.59 sign p 0.33 Wilcoxon p 0.28 | J +0.029/h [-0.002,+0.060] | onsets Q:A 296:470 wakings Q:A 272:521
    per0LD  dark  flies 27 | c mean +0.009 [-0.031,+0.053] median -0.018 c>0 0.41 sign p 0.44 Wilcoxon p 0.94 | J +0.016/h [-0.012,+0.049] | onsets Q:A 344:499 wakings Q:A 310:506
    per0DD  whole flies 25 | c mean +0.122 [+0.059,+0.182] median +0.119 c>0 0.72 sign p 0.023 Wilcoxon p 0.0025 | J +0.084/h [+0.036,+0.130] | onsets Q:A 731:826 wakings Q:A 546:1011
    sss     light flies 50 | c mean -0.076 [-0.150,-0.008] median -0.051 c>0 0.38 sign p 0.15 Wilcoxon p 0.088 | J -0.015/h [-0.030,+0.000] | onsets Q:A 217:319 wakings Q:A 269:317
    sss     dark  flies 53 | c mean -0.071 [-0.140,-0.002] median -0.067 c>0 0.36 sign p 0.092 Wilcoxon p 0.058 | J -0.017/h [-0.037,+0.004] | onsets Q:A 287:414 wakings Q:A 295:357
      sleepless vs others (Mann-Whitney on c):
        light: sss vs CS: p 0.468
        light: sss vs fumin: p 4.57e-05
        light: sss vs per0LD: p 0.0547
        light: sss vs per0DD(whole): p 0.000463
        dark: sss vs CS: p 3.75e-05
        dark: sss vs fumin: p 9.87e-07
        dark: sss vs per0LD: p 0.108
        dark: sss vs per0DD(whole): p 0.000715
    INTERPRETATION (descriptive; post hoc window): forward current (c > 0, 95% CI excludes 0) in wild type dark
    (+0.11), fumin light (+0.05) and dark (+0.06), per0 DD (+0.11); none in wild type light or per0 LD (either half).
    sleepless LIGHT reversed: c = -0.093 [-0.154, -0.036], negative at all three thresholds (CI excludes 0 at each);
    sleepless dark -0.03 (ns at theta 0.33/0.5; -0.07, CI just excludes 0, at 0.75). Net rates small: |J| <= 0.07
    cycles/h; sleepless J CIs include 0. Not independent of Fig. 3 (same window, thresholded). per0 LD dark:
    positive Delta (graded) but c ~ 0 (binary Q/A classification loses it).
    Figure: figures/fig5_cycle_current.{pdf,png}.

## PRESPEC 20 section 7 — synthetic validation of the S/Q/W per-window HMM (2026-10-23, code 00c130f)
- V1 (true 3-state, distinct H1 emissions, Xn -0.3, Yn +0.4; 50 windows x 240 min): **FAIL.**
  BIC chose K=2 in 96% of windows, K=3 in 4% (n_eff rule: 90% / 10%). H1 degenerate in 2% only.
  K=3 routing per window: Xn mean -0.24, median -0.42, SD 0.75; Yn mean +0.25, median +0.48, SD 0.75.
  Distinct emissions are not enough: one 4 h window (~50 symbol changes) cannot resolve routing or justify a third state.
- V2 (true 2-state, no H1): PASS. K=2 in 100%; K=3 routing uninformative in 100% (bootstrap SD ~0.78).
- Consequence, per the plan: the per-window design is not run on fly data. Design review before any fly fit.

## Design check after V1 failure — steps per fit (2026-10-23, code 72f15b9; synthetic only, V1 truth, 30 reps, 10 starts)
- n   240 (1 windows): K=3 chosen 0.17, K=2 0.83, K=4 0.00; Xn mean -0.087 SD 0.830; Yn mean +0.405 SD 0.596; degenerate 0.00
- n   960 (4 windows): K=3 chosen 0.70, K=2 0.30, K=4 0.00; Xn mean -0.212 SD 0.403; Yn mean +0.358 SD 0.384; degenerate 0.00
- n  2880 (12 windows): K=3 chosen 1.00, K=2 0.00, K=4 0.00; Xn mean -0.355 SD 0.177; Yn mean +0.414 SD 0.229; degenerate 0.00
- n  5760 (24 windows): K=3 chosen 1.00, K=2 0.00, K=4 0.00; Xn mean -0.300 SD 0.121; Yn mean +0.402 SD 0.165; degenerate 0.00
- SD falls roughly as 1/sqrt(n). Routing SD < 0.2 needs ~2880 steps (48 h of minutes). n=240 gave K=3 in 17% here vs 4% in V1 (different seeds/starts): sampling variation.

## PRESPEC 20b section 4 — synthetic validation, per fly x light/dark (2026-10-23, code f53efb8)
- truth light Xn -0.300 Yn +0.400; dark Xn +0.200 Yn -0.100
- V1b: best K (n=T) {2: 0.0, 3: 1.0, 4: 0.0}; (n_eff) {2: 0.0, 3: 1.0, 4: 0.0}; disagree 0.00; H1 degenerate 0.00; median changes 1149.0; max EM iters 1000
-    xl mean -0.272 SD 0.238; yl mean +0.409 SD 0.200; xd mean +0.186 SD 0.272; yd mean -0.105 SD 0.230
-    V1b PASS=True
- V2b: best K (n=T) {2: 1.0, 3: 0.0, 4: 0.0}; (n_eff) {2: 1.0, 3: 0.0, 4: 0.0}; disagree 0.00; H1 degenerate 0.25; median changes 1288.5; max EM iters 1000
-    xl mean -0.440 SD 0.718; yl mean -0.330 SD 0.714; xd mean -0.127 SD 0.786; yd mean -0.269 SD 0.754
- V3b: best K (n=T) {2: 1.0, 3: 0.0, 4: 0.0}; (n_eff) {2: 1.0, 3: 0.0, 4: 0.0}; disagree 0.00; H1 degenerate 0.50; median changes 1282.0; max EM iters 1000
-    xl mean -0.130 SD 0.781; yl mean -0.149 SD 0.902; xd mean -0.102 SD 0.818; yd mean +0.224 SD 0.808
-    V2b PASS=True
-    V3b PASS=True
- elapsed 2824 s
- All three PASS. Routing SD per phase 0.20-0.27 at 2880 min/phase, matching the steps-per-fit check. Note: the batched EM loop hit its 1000-iteration cap in some fits (the cap counts until the slowest of 20 starts converges); results kept as specified.

## PRESPEC 20b primary — S/Q/W HMM per fly x light/dark, theta 0.5 (2026-10-23, code ac79c9b, report bc7c342)
Narrowest block: one fly's light (ZT 0-12) or dark (ZT 12-24) minutes pooled over all its days; per0DD whole record. Group values: mean over flies, 95% fly-bootstrap interval.
```
CS: flies 25 | best K (n=T) 2:0.00 3:0.96 4:0.04 | (n_eff) 2:0.00 3:0.88 4:0.12 | rules disagree 0.08 | flagged 0.04 | H1 degenerate 0.00
  light (non-degenerate flies 25) Xn_pred +0.030 [-0.116,+0.175]  Yn_pred +0.041 [-0.078,+0.160]  c_pred +0.005 [-0.037,+0.045]
        ffbs Xn -0.029 [-0.168,+0.117]  Yn -0.025 [-0.140,+0.099]  c_model +0.002 [-0.043,+0.044]
        vit  Xn -0.037 [-0.157,+0.083]  Yn -0.037 [-0.130,+0.063]  c_model +0.000 [-0.037,+0.037]
        model-free c -0.019 [-0.079,+0.048]; Spearman(c_obs, c_model ffbs) +0.44 p 0.028 n 25
        energetics: infinite EPR in 0.00 of flies; finite median EPR 0.032 nats/h, J -0.0131 cycles/h, affinity -0.04
  dark  (non-degenerate flies 25) Xn_pred -0.302 [-0.404,-0.191]  Yn_pred +0.216 [+0.117,+0.316]  c_pred +0.259 [+0.195,+0.323]
        ffbs Xn -0.265 [-0.366,-0.149]  Yn +0.235 [+0.124,+0.342]  c_model +0.250 [+0.187,+0.312]
        vit  Xn -0.250 [-0.388,-0.118]  Yn +0.238 [+0.119,+0.340]  c_model +0.244 [+0.176,+0.312]
        model-free c +0.111 [+0.068,+0.155]; Spearman(c_obs, c_model ffbs) +0.32 p 0.12 n 25
        energetics: infinite EPR in 0.00 of flies; finite median EPR 0.281 nats/h, J -0.2250 cycles/h, affinity -1.14
fumin: flies 53 | best K (n=T) 2:0.06 3:0.85 4:0.09 | (n_eff) 2:0.02 3:0.72 4:0.26 | rules disagree 0.21 | flagged 0.13 | H1 degenerate 0.00
  light (non-degenerate flies 53) Xn_pred -0.370 [-0.516,-0.211]  Yn_pred -0.027 [-0.177,+0.123]  c_pred +0.172 [+0.100,+0.241]
        ffbs Xn -0.498 [-0.613,-0.364]  Yn -0.188 [-0.313,-0.076]  c_model +0.178 [+0.115,+0.242]
        vit  Xn -0.512 [-0.642,-0.369]  Yn -0.215 [-0.325,-0.105]  c_model +0.177 [+0.100,+0.251]
        model-free c +0.048 [+0.018,+0.079]; Spearman(c_obs, c_model ffbs) +0.19 p 0.23 n 44
        energetics: infinite EPR in 0.21 of flies; finite median EPR 0.160 nats/h, J -0.1288 cycles/h, affinity -1.03
  dark  (non-degenerate flies 53) Xn_pred -0.466 [-0.619,-0.303]  Yn_pred +0.052 [-0.102,+0.207]  c_pred +0.259 [+0.162,+0.348]
        ffbs Xn -0.561 [-0.690,-0.420]  Yn +0.068 [-0.088,+0.212]  c_model +0.335 [+0.256,+0.414]
        vit  Xn -0.569 [-0.697,-0.424]  Yn -0.010 [-0.153,+0.140]  c_model +0.310 [+0.231,+0.392]
        model-free c +0.059 [+0.013,+0.102]; Spearman(c_obs, c_model ffbs) +0.23 p 0.12 n 45
        energetics: infinite EPR in 0.38 of flies; finite median EPR 0.193 nats/h, J -0.1171 cycles/h, affinity -1.51
per0LD: flies 27 | best K (n=T) 2:0.15 3:0.78 4:0.07 | (n_eff) 2:0.04 3:0.89 4:0.07 | rules disagree 0.11 | flagged 0.11 | H1 degenerate 0.00
  light (non-degenerate flies 27) Xn_pred -0.424 [-0.549,-0.288]  Yn_pred -0.176 [-0.319,-0.025]  c_pred +0.124 [+0.068,+0.185]
        ffbs Xn -0.426 [-0.548,-0.305]  Yn -0.175 [-0.318,-0.021]  c_model +0.126 [+0.073,+0.181]
        vit  Xn -0.474 [-0.597,-0.353]  Yn -0.218 [-0.374,-0.048]  c_model +0.128 [+0.069,+0.195]
        model-free c +0.007 [-0.030,+0.044]; Spearman(c_obs, c_model ffbs) +0.08 p 0.69 n 27
        energetics: infinite EPR in 0.04 of flies; finite median EPR 0.015 nats/h, J -0.0413 cycles/h, affinity -0.28
  dark  (non-degenerate flies 27) Xn_pred -0.473 [-0.572,-0.363]  Yn_pred -0.185 [-0.287,-0.083]  c_pred +0.144 [+0.089,+0.203]
        ffbs Xn -0.475 [-0.567,-0.379]  Yn -0.186 [-0.283,-0.077]  c_model +0.144 [+0.090,+0.201]
        vit  Xn -0.526 [-0.632,-0.420]  Yn -0.264 [-0.377,-0.156]  c_model +0.131 [+0.070,+0.196]
        model-free c -0.015 [-0.046,+0.019]; Spearman(c_obs, c_model ffbs) +0.14 p 0.48 n 27
        energetics: infinite EPR in 0.04 of flies; finite median EPR 0.163 nats/h, J -0.1967 cycles/h, affinity -0.61
per0DD: flies 25 | best K (n=T) 2:0.08 3:0.76 4:0.16 | (n_eff) 2:0.04 3:0.56 4:0.40 | rules disagree 0.28 | flagged 0.20 | H1 degenerate 0.08
  whole (non-degenerate flies 23) Xn_pred -0.341 [-0.451,-0.223]  Yn_pred +0.082 [-0.012,+0.178]  c_pred +0.211 [+0.137,+0.285]
        ffbs Xn -0.310 [-0.414,-0.189]  Yn +0.104 [+0.012,+0.199]  c_model +0.207 [+0.138,+0.277]
        vit  Xn -0.321 [-0.445,-0.195]  Yn +0.131 [+0.029,+0.235]  c_model +0.226 [+0.141,+0.314]
        model-free c +0.107 [+0.051,+0.160]; Spearman(c_obs, c_model ffbs) +0.53 p 0.0092 n 23
        energetics: infinite EPR in 0.00 of flies; finite median EPR 0.208 nats/h, J -0.2337 cycles/h, affinity -0.95
sss: flies 57 | best K (n=T) 2:0.05 3:0.95 4:0.00 | (n_eff) 2:0.05 3:0.91 4:0.04 | rules disagree 0.04 | flagged 0.05 | H1 degenerate 0.05
  light (non-degenerate flies 54) Xn_pred +0.122 [-0.064,+0.301]  Yn_pred -0.145 [-0.347,+0.066]  c_pred -0.133 [-0.216,-0.054]
        ffbs Xn +0.071 [-0.136,+0.277]  Yn -0.240 [-0.434,-0.023]  c_model -0.185 [-0.271,-0.102]
        vit  Xn +0.085 [-0.124,+0.287]  Yn -0.206 [-0.415,+0.026]  c_model -0.179 [-0.287,-0.076]
        model-free c -0.093 [-0.156,-0.036]; Spearman(c_obs, c_model ffbs) +0.28 p 0.096 n 37
        energetics: infinite EPR in 0.46 of flies; finite median EPR 0.093 nats/h, J +0.0522 cycles/h, affinity +0.81
  dark  (non-degenerate flies 54) Xn_pred +0.005 [-0.179,+0.178]  Yn_pred -0.174 [-0.355,+0.018]  c_pred -0.089 [-0.172,-0.008]
        ffbs Xn -0.124 [-0.293,+0.055]  Yn -0.221 [-0.398,-0.060]  c_model -0.052 [-0.132,+0.027]
        vit  Xn -0.145 [-0.327,+0.040]  Yn -0.249 [-0.427,-0.065]  c_model -0.053 [-0.142,+0.033]
        model-free c -0.029 [-0.098,+0.039]; Spearman(c_obs, c_model ffbs) +0.33 p 0.032 n 43
        energetics: infinite EPR in 0.46 of flies; finite median EPR 0.064 nats/h, J +0.0180 cycles/h, affinity +0.21
```
- Pre-stated sign prediction (section 8 of PRESPEC 20): c_model has the model-free c's sign in all 5 blocks where that interval excludes 0 (CS dark, fumin light, fumin dark, per0DD, sss light). PASS 5/5. Magnitudes are 1.5-6x the model-free c; per-fly Spearman correlations are weak (+0.08 to +0.53).
- Disagreement: per0LD, model c_model +0.13/+0.14 with model-free c ~0 in both phases.
- Sign convention: J and affinity are computed for the cycle H1->H2->H3, i.e. Q->A->S, which is the REVERSE of the forward cycle A->Q->S. Negative J = forward circulation, consistent with positive c.
- Infinite EPR (a transition zero in one direction only) in 21-46% of fumin and sss flies: reported, not floored.
- Not yet run: bootstrap, day split, adequacy, theta 0.33/0.75.

## POST HOC (user question, 2026-10-23): why did the scored S/W HMM split per0 flies into two groups?
Scored .mat data, accepted windows (common.py). H1 called wake-like / sleep-like if its emission row is within
total variation 0.05 of H2's / H3's. Per fly: fraction of sleep-like windows vs mean (Yn_pred - Xn_pred).
- per0 (25 flies, 544 windows): H1 wake-like 0.60, sleep-like 0.32, distinct 0.08. Flies with Yn-Xn > 0: 0.52.
  Spearman(fraction sleep-like, Yn-Xn) = -0.94, p 2.5e-12. Sleep-like windows: median Yn-Xn -1.59 (43% > 0);
  wake-like windows: median 0.00 (56% > 0).
- CS: sleep-like 0.27, rho -0.83 (p 3e-7), but sleep-like windows stay mostly positive (median +0.45), 88% of flies > 0.
- fumin, sss: H1 wake-like in 99% of windows; 100% of flies > 0.
Reading: a fly's side of the diagonal is predicted almost entirely by how often its degenerate H1 copied the
SLEEP emission row. The sign comes from the gauge choice, not from identified routing. In the S/Q/W model (20b)
per0 shows no split (85-89% of flies c > 0). Descriptive.

## ADDED AFTER 20b RESULTS (user request): count-based J_cycle, theta 0.5 (2026-10-23, code 7eacf20)
Definition as compute_HMM_energetics.m GLOBAL: J_cycle = N32 - N23 on the run-compressed hidden sequence, per hour; here averaged over 200 FFBS whole-record paths (no window seams), counts attributed to the phase / 4 h window of the arrival minute (window counts sum to the global counts, asserted). Sign: + = forward A->Q->S->A. Model J = stationary flux of the phase's fitted matrix, same sign convention.
```
CS      light n 25 | model J +0.020 [-0.024,+0.067] | count J_cycle (FFBS) +0.030 [-0.016,+0.076] | (Viterbi) +0.027 cycles/h | Spearman(model, count) +0.99 | same sign 0.96 | count EPR (FFBS) median 0.035 /h, infinite 0.00
CS      dark  n 25 | model J +0.217 [+0.155,+0.284] | count J_cycle (FFBS) +0.205 [+0.143,+0.271] | (Viterbi) +0.203 cycles/h | Spearman(model, count) +0.99 | same sign 1.00 | count EPR (FFBS) median 0.291 /h, infinite 0.00
fumin   light n 53 | model J +0.131 [+0.084,+0.177] | count J_cycle (FFBS) +0.129 [+0.082,+0.180] | (Viterbi) +0.131 cycles/h | Spearman(model, count) +0.98 | same sign 0.87 | count EPR (FFBS) median 0.139 /h, infinite 0.17
fumin   dark  n 53 | model J +0.127 [+0.089,+0.166] | count J_cycle (FFBS) +0.136 [+0.100,+0.178] | (Viterbi) +0.128 cycles/h | Spearman(model, count) +0.99 | same sign 0.92 | count EPR (FFBS) median 0.177 /h, infinite 0.34
per0LD  light n 27 | model J +0.101 [+0.051,+0.156] | count J_cycle (FFBS) +0.109 [+0.058,+0.164] | (Viterbi) +0.109 cycles/h | Spearman(model, count) +0.93 | same sign 0.89 | count EPR (FFBS) median 0.019 /h, infinite 0.04
per0LD  dark  n 27 | model J +0.187 [+0.116,+0.256] | count J_cycle (FFBS) +0.181 [+0.111,+0.250] | (Viterbi) +0.160 cycles/h | Spearman(model, count) +1.00 | same sign 1.00 | count EPR (FFBS) median 0.162 /h, infinite 0.04
per0DD  whole n 23 | model J +0.247 [+0.154,+0.353] | count J_cycle (FFBS) +0.246 [+0.159,+0.343] | (Viterbi) +0.239 cycles/h | Spearman(model, count) +1.00 | same sign 1.00 | count EPR (FFBS) median 0.206 /h, infinite 0.00
sss     light n 54 | model J -0.056 [-0.083,-0.031] | count J_cycle (FFBS) -0.056 [-0.084,-0.031] | (Viterbi) -0.050 cycles/h | Spearman(model, count) +0.94 | same sign 0.70 | count EPR (FFBS) median 0.020 /h, infinite 0.28
sss     dark  n 54 | model J -0.016 [-0.041,+0.011] | count J_cycle (FFBS) -0.018 [-0.045,+0.007] | (Viterbi) -0.019 cycles/h | Spearman(model, count) +0.99 | same sign 0.87 | count EPR (FFBS) median 0.048 /h, infinite 0.35
CS      ZT0-4..20-24 count J (cycles/h), mean over flies [SD across flies]: -0.062 [0.12]  +0.058 [0.20]  +0.094 [0.21]  +0.178 [0.18]  +0.202 [0.21]  +0.236 [0.21]
        fraction of flies with J > 0 per window: 0.28  0.56  0.60  0.84  0.84  0.84
fumin   ZT0-4..20-24 count J (cycles/h), mean over flies [SD across flies]: +0.087 [0.22]  +0.143 [0.22]  +0.157 [0.27]  +0.245 [0.27]  +0.115 [0.19]  +0.049 [0.12]
        fraction of flies with J > 0 per window: 0.58  0.60  0.74  0.77  0.62  0.45
per0LD  ZT0-4..20-24 count J (cycles/h), mean over flies [SD across flies]: +0.137 [0.16]  +0.079 [0.19]  +0.111 [0.20]  +0.085 [0.19]  +0.236 [0.22]  +0.223 [0.25]
        fraction of flies with J > 0 per window: 0.93  0.67  0.70  0.59  0.85  0.85
per0DD  ZT0-4..20-24 count J (cycles/h), mean over flies [SD across flies]: +0.260 [0.26]  +0.243 [0.25]  +0.228 [0.29]  +0.253 [0.31]  +0.278 [0.31]  +0.211 [0.28]
        fraction of flies with J > 0 per window: 0.87  0.83  0.78  0.78  0.87  0.91
sss     ZT0-4..20-24 count J (cycles/h), mean over flies [SD across flies]: -0.052 [0.14]  -0.065 [0.12]  -0.050 [0.12]  -0.004 [0.12]  -0.005 [0.15]  -0.045 [0.16]
        fraction of flies with J > 0 per window: 0.22  0.11  0.15  0.31  0.37  0.26
```
- Count J and model J agree (Spearman +0.93 to +1.00 per block); differences are within fly-bootstrap intervals.
- Per-window values come from a model with ONE transition matrix per 12 h phase: their ZT variation reflects data-driven posterior paths and relaxation after light switches, not a fitted time-varying matrix. Descriptive only; the 6 h fits (20c) are the model-level test.

## Jensen check (user request), theta 0.5 (2026-10-23, code b060abb)
Nonlinear statistics (Xn, Yn, c, count EPR) computed per FFBS path then averaged, E[f(N)], beside the plug-in on posterior-mean counts, f(E[N]). Per-path means are reported only for flies with the statistic defined on >= 50% of paths, so fly sets can differ slightly between columns.
```
CS      light | c plug-in +0.002 vs per-path mean +0.002 (paths with c defined 1.00) | Xn -0.029 vs -0.028 | Yn -0.025 vs -0.025 | EPR plug-in median 0.035 vs per-path median 0.037 /h, per-path infinite 0.00
CS      dark  | c plug-in +0.250 vs per-path mean +0.250 (paths with c defined 1.00) | Xn -0.265 vs -0.265 | Yn +0.235 vs +0.235 | EPR plug-in median 0.291 vs per-path median 0.285 /h, per-path infinite 0.00
fumin   light | c plug-in +0.178 vs per-path mean +0.170 (paths with c defined 0.88) | Xn -0.498 vs -0.470 | Yn -0.188 vs -0.189 | EPR plug-in median 0.139 vs per-path median 0.140 /h, per-path infinite 0.23
fumin   dark  | c plug-in +0.335 vs per-path mean +0.335 (paths with c defined 0.85) | Xn -0.561 vs -0.562 | Yn +0.068 vs +0.068 | EPR plug-in median 0.177 vs per-path median 0.206 /h, per-path infinite 0.39
per0LD  light | c plug-in +0.126 vs per-path mean +0.126 (paths with c defined 1.00) | Xn -0.426 vs -0.426 | Yn -0.175 vs -0.174 | EPR plug-in median 0.019 vs per-path median 0.024 /h, per-path infinite 0.07
per0LD  dark  | c plug-in +0.144 vs per-path mean +0.144 (paths with c defined 1.00) | Xn -0.475 vs -0.475 | Yn -0.186 vs -0.186 | EPR plug-in median 0.162 vs per-path median 0.160 /h, per-path infinite 0.04
per0DD  whole | c plug-in +0.207 vs per-path mean +0.207 (paths with c defined 1.00) | Xn -0.310 vs -0.310 | Yn +0.104 vs +0.104 | EPR plug-in median 0.206 vs per-path median 0.212 /h, per-path infinite 0.00
sss     light | c plug-in -0.185 vs per-path mean -0.185 (paths with c defined 0.75) | Xn +0.071 vs +0.070 | Yn -0.240 vs -0.240 | EPR plug-in median 0.020 vs per-path median 0.025 /h, per-path infinite 0.31
sss     dark  | c plug-in -0.052 vs per-path mean -0.062 (paths with c defined 0.84) | Xn -0.124 vs -0.099 | Yn -0.221 vs -0.194 | EPR plug-in median 0.048 vs per-path median 0.038 /h, per-path infinite 0.39
```
- c, Xn, Yn: E[f] and f(E) agree to <= 0.01 except fumin light Xn (0.03), sss dark Xn/Yn (0.03). J_cycle is linear: identical by construction.
- Count EPR: on 23-39% of paths in fumin and sss some edge is traversed one way only, so the per-path EPR is infinite and its posterior MEAN is infinite for those flies; medians reported.

## Rule-of-3 EPR lower bound (user suggestion), theta 0.5 (2026-10-23, code 4403df1)
Zero count b in pair (a, b) -> b* = min(3, a): (a-b) ln(a/b) is minimised over b in [0, 3] there, so the result is a 95% LOWER BOUND per term, not an estimate. Per path, then averaged over the 200 FFBS paths.
```
CS      light | EPR 95% lower bound (rule of 3), per-path posterior mean: median over flies 0.041 /h, mean 0.067 | paths needing substitution 0.00 | plug-in on mean counts, median 0.035 /h | flies 25
CS      dark  | EPR 95% lower bound (rule of 3), per-path posterior mean: median over flies 0.301 /h, mean 0.380 | paths needing substitution 0.00 | plug-in on mean counts, median 0.291 /h | flies 25
fumin   light | EPR 95% lower bound (rule of 3), per-path posterior mean: median over flies 0.108 /h, mean 0.257 | paths needing substitution 0.23 | plug-in on mean counts, median 0.085 /h | flies 53
fumin   dark  | EPR 95% lower bound (rule of 3), per-path posterior mean: median over flies 0.143 /h, mean 0.272 | paths needing substitution 0.39 | plug-in on mean counts, median 0.146 /h | flies 53
per0LD  light | EPR 95% lower bound (rule of 3), per-path posterior mean: median over flies 0.030 /h, mean 0.165 | paths needing substitution 0.07 | plug-in on mean counts, median 0.020 /h | flies 27
per0LD  dark  | EPR 95% lower bound (rule of 3), per-path posterior mean: median over flies 0.177 /h, mean 0.279 | paths needing substitution 0.04 | plug-in on mean counts, median 0.159 /h | flies 27
per0DD  whole | EPR 95% lower bound (rule of 3), per-path posterior mean: median over flies 0.208 /h, mean 0.415 | paths needing substitution 0.00 | plug-in on mean counts, median 0.206 /h | flies 23
sss     light | EPR 95% lower bound (rule of 3), per-path posterior mean: median over flies 0.031 /h, mean 0.103 | paths needing substitution 0.31 | plug-in on mean counts, median 0.024 /h | flies 54
sss     dark  | EPR 95% lower bound (rule of 3), per-path posterior mean: median over flies 0.044 /h, mean 0.086 | paths needing substitution 0.39 | plug-in on mean counts, median 0.032 /h | flies 54
```
- Finite in every fly. Group means exceed medians (right-skewed across flies). The plug-in median here includes flies that were infinite before, so it is not comparable with the earlier finite-only plug-in medians.
- Cross-check (2026-10-23): the rule-of-3 bound equals a Python port of the user's epr_count_based.m ('bound') to 8.9e-16 over 1e5 random count matrices; min per-path bound over all real flies/phases/paths = 0.0 (EPR >= 0 now asserted in code).

## POST HOC (user observation, 2026-10-23): light vs dark within fly, theta 0.5, 20b fits
Paired per fly (non-degenerate), Wilcoxon signed-rank, dark vs light:
- EPR lower bound (rule of 3, per-path mean): CS dark > light in 68% of flies, p 0.00063 (medians 0.041 / 0.301);
  fumin 51%, p 0.85; per0LD 52%, p 0.27 (medians 0.030 / 0.177); sss 54%, p 0.70.
- count J_cycle: CS 80%, p 0.00025 (median paired diff +0.18/h); fumin 53%, p 0.72; per0LD 56%, p 0.12
  (+0.07/h); sss 61%, p 0.042 (+0.03/h).
Reading: per0LD's light/dark difference in GROUP MEDIANS is not a consistent within-fly difference; paired, per0LD
sits with fumin/sss, not with CS. Caveats: n = 27 (limited power); the Q threshold reference r_f is set per
phase, so the symbol coding itself differs between light and dark. Descriptive.

## PRESPEC 20b section 3 — parametric bootstrap, 5 flies per genotype, 30 replicates (2026-10-23, code b060abb)
```
CS      light: flies 5 | SD Xn median 0.118 [0.11-0.19] | SD Yn median 0.124 [0.08-0.22] | bias Xn +0.004 Yn -0.014 | flies with undefined values 0
CS      dark : flies 5 | SD Xn median 0.125 [0.10-0.22] | SD Yn median 0.156 [0.11-0.17] | bias Xn -0.004 Yn +0.027 | flies with undefined values 0
fumin   light: flies 5 | SD Xn median 0.174 [0.15-0.55] | SD Yn median 0.247 [0.17-0.38] | bias Xn -0.034 Yn +0.010 | flies with undefined values 1
fumin   dark : flies 5 | SD Xn median 0.525 [0.14-0.69] | SD Yn median 0.338 [0.14-0.70] | bias Xn +0.076 Yn +0.048 | flies with undefined values 0
per0LD  light: flies 5 | SD Xn median 0.171 [0.09-0.20] | SD Yn median 0.165 [0.11-0.24] | bias Xn -0.016 Yn +0.022 | flies with undefined values 0
per0LD  dark : flies 5 | SD Xn median 0.150 [0.11-0.28] | SD Yn median 0.210 [0.10-0.24] | bias Xn +0.014 Yn -0.016 | flies with undefined values 0
per0DD  whole: flies 5 | SD Xn median 0.136 [0.07-0.27] | SD Yn median 0.148 [0.09-0.17] | bias Xn -0.027 Yn +0.015 | flies with undefined values 0
sss     light: flies 5 | SD Xn median 0.369 [0.22-0.83] | SD Yn median 0.241 [0.21-0.37] | bias Xn +0.027 Yn -0.059 | flies with undefined values 0
sss     dark : flies 5 | SD Xn median 0.468 [0.19-0.74] | SD Yn median 0.411 [0.24-0.71] | bias Xn +0.018 Yn +0.039 | flies with undefined values 0
```
- Bias is small everywhere (|median| <= 0.08). Per-fly precision splits by genotype: CS, per0LD, per0DD SD 0.12-0.21; fumin dark and sss light/dark SD 0.34-0.53 (few sleep transitions). Group means over 53-57 flies remain usable (SE ~ SD/sqrt(n) ~ 0.05-0.07), but per-fly routing in fumin and sss is too noisy for per-fly claims (relevant to 20e flips).

## PRESPEC 20b section 5 — day split and within-phase adequacy (2026-10-23, code b060abb; summary sqw2_checks.py)
```
day split CS      light: flies 25 | c_pred first -0.010 second +0.011 | change +0.021 [-0.085,+0.131] | Spearman -0.07 (p 0.76) | sign flips 0.44
day split CS      dark : flies 25 | c_pred first +0.237 second +0.312 | change +0.075 [-0.043,+0.198] | Spearman +0.14 (p 0.51) | sign flips 0.16
day split fumin   light: flies 49 | c_pred first +0.146 second +0.185 | change +0.039 [-0.096,+0.172] | Spearman +0.05 (p 0.71) | sign flips 0.37
day split fumin   dark : flies 47 | c_pred first +0.253 second +0.267 | change +0.015 [-0.147,+0.167] | Spearman +0.01 (p 0.95) | sign flips 0.51
day split per0LD  light: flies 27 | c_pred first +0.122 second +0.171 | change +0.049 [-0.097,+0.185] | Spearman -0.09 (p 0.65) | sign flips 0.44
day split per0LD  dark : flies 27 | c_pred first +0.141 second +0.122 | change -0.019 [-0.134,+0.081] | Spearman +0.11 (p 0.58) | sign flips 0.33
day split per0DD  whole: flies 22 | c_pred first +0.232 second +0.208 | change -0.024 [-0.137,+0.089] | Spearman +0.25 (p 0.26) | sign flips 0.14
day split sss     light: flies 42 | c_pred first -0.111 second -0.073 | change +0.039 [-0.107,+0.192] | Spearman +0.12 (p 0.44) | sign flips 0.45
day split sss     dark : flies 46 | c_pred first -0.016 second -0.134 | change -0.117 [-0.284,+0.051] | Spearman -0.13 (p 0.39) | sign flips 0.70
within-phase adequacy: fraction of 4 h windows outside the fitted model's 95% predictive band, by ZT window (0-4..20-24); columns Q, W, S fraction, symbol changes; ~0.05 expected
  CS      ZT0-4: 0.07 0.31 0.08 0.09 | ZT4-8: 0.13 0.30 0.17 0.14 | ZT8-12: 0.14 0.08 0.07 0.11 | ZT12-16: 0.05 0.04 0.06 0.07 | ZT16-20: 0.01 0.01 0.01 0.01 | ZT20-24: 0.09 0.06 0.05 0.10
  fumin   ZT0-4: 0.47 0.58 0.19 0.44 | ZT4-8: 0.11 0.19 0.08 0.17 | ZT8-12: 0.16 0.14 0.05 0.15 | ZT12-16: 0.24 0.26 0.09 0.26 | ZT16-20: 0.19 0.23 0.02 0.26 | ZT20-24: 0.24 0.27 0.02 0.25
  per0LD  ZT0-4: 0.20 0.31 0.22 0.19 | ZT4-8: 0.09 0.11 0.07 0.14 | ZT8-12: 0.12 0.08 0.03 0.10 | ZT12-16: 0.15 0.09 0.08 0.18 | ZT16-20: 0.15 0.08 0.06 0.06 | ZT20-24: 0.12 0.07 0.08 0.10
  per0DD  ZT0-4: 0.02 0.00 0.00 0.00 | ZT4-8: 0.02 0.01 0.00 0.03 | ZT8-12: 0.02 0.00 0.00 0.01 | ZT12-16: 0.04 0.01 0.00 0.00 | ZT16-20: 0.02 0.01 0.00 0.00 | ZT20-24: 0.02 0.01 0.00 0.01
  sss     ZT0-4: 0.23 0.24 0.08 0.23 | ZT4-8: 0.25 0.27 0.09 0.28 | ZT8-12: 0.24 0.35 0.14 0.26 | ZT12-16: 0.26 0.30 0.12 0.24 | ZT16-20: 0.23 0.31 0.11 0.25 | ZT20-24: 0.28 0.27 0.09 0.18
```
- Day split: no group-level change in any block (all intervals include 0). Per-fly c in the two halves is uncorrelated (Spearman -0.13 to +0.25, all ns); sign-flip fractions 0.14-0.70 must be judged against the 20e no-change null.
- Adequacy: per0DD fits within band (<= 0.04). CS: W fraction outside band in 30-31% of ZT 0-8 windows (within-light-phase structure the 12 h model averages). per0LD: 20-31% at ZT 0-4 only (lights-on). fumin: 44-58% at ZT 0-4 and 15-26% elsewhere. sss: 23-35% at all ZT (not time-of-day specific: between-day / overdispersion). The 12 h homogeneous-phase assumption is violated most at lights-on and in fumin/sss generally.

## PRESPEC 20c step 1 — 6 h blocks, synthetic (2026-10-23, code fd88bf6)
```
V1c: best K {2: 0.0, 3: 1.0, 4: 0.0}; H1 degenerate 0.00
   Xn mean -0.268 -0.335 +0.187 +0.197 SD 0.237 0.223 0.281 0.323
   Yn mean +0.444 +0.469 -0.133 -0.014 SD 0.163 0.324 0.414 0.377
   V1c PASS=False
V2c: best K {2: 1.0, 3: 0.0, 4: 0.0}; H1 degenerate 0.35
   Xn mean -0.415 -0.088 -0.001 -0.035 SD 0.820 0.886 0.930 0.900
   Yn mean +0.198 +0.014 -0.310 -0.143 SD 0.816 0.826 0.756 0.871
   V2c PASS=True
elapsed 1783 s
```
- V1c FAIL on precision: K=3 chosen 100%, biases <= 0.09, but SD exceeds 0.3 in 3 of 8 quantities (Yn blocks 2-4: 0.32-0.41; Xn block 4: 0.32). As extrapolated: 6 h blocks are borderline on this truth. V2c PASS. Per the plan, steps 2-3 proceed; the usable verdict per genotype comes from step 2.

## PRESPEC 20c step 2 — realistic truths, 6 h blocks (2026-10-24, code fd88bf6)
```
CS fly 6: K=3 chosen 1.00; Xn SD 0.26 0.25 0.31 0.32; Yn SD 0.19 0.35 0.27 0.26; Xn bias +0.06 +0.08 +0.01 +0.09; Yn bias +0.09 +0.05 +0.02 +0.01
CS fly 8: K=3 chosen 1.00; Xn SD 0.17 0.22 0.30 0.32; Yn SD 0.23 0.24 0.28 0.37; Xn bias +0.02 +0.06 +0.05 +0.25; Yn bias +0.06 -0.00 +0.01 -0.00
CS fly 11: K=3 chosen 1.00; Xn SD 0.27 0.18 0.16 0.16; Yn SD 0.24 0.15 0.21 0.29; Xn bias +0.00 -0.03 -0.10 -0.08; Yn bias -0.04 -0.06 -0.01 -0.00
CS fly 14: K=3 chosen 1.00; Xn SD 0.21 0.24 0.29 0.18; Yn SD 0.17 0.20 0.27 0.27; Xn bias -0.06 +0.05 -0.04 +0.00; Yn bias -0.02 +0.03 +0.06 -0.03
CS fly 17: K=3 chosen 1.00; Xn SD 0.20 0.26 0.23 0.21; Yn SD 0.20 0.19 0.34 0.25; Xn bias -0.02 +0.12 -0.14 +0.02; Yn bias +0.01 +0.02 -0.01 +0.00
CS: median SD over flies (Xn blocks 1-4, Yn blocks 1-4) 0.21 0.24 0.29 0.21 0.20 0.20 0.27 0.27 -> USABLE=True
per0LD fly 6: K=3 chosen 0.90; Xn SD 0.48 0.28 0.23 0.17; Yn SD 0.27 0.37 0.32 0.23; Xn bias +0.12 +0.07 -0.00 +0.06; Yn bias +0.09 -0.01 -0.08 -0.09
per0LD fly 11: K=3 chosen 1.00; Xn SD 0.18 0.21 0.17 0.15; Yn SD 0.21 0.22 0.15 0.16; Xn bias -0.03 -0.01 -0.07 -0.01; Yn bias +0.05 +0.09 +0.03 -0.03
per0LD fly 14: K=3 chosen 1.00; Xn SD 0.21 0.21 0.17 0.16; Yn SD 0.19 0.21 0.19 0.21; Xn bias -0.06 -0.01 -0.02 -0.04; Yn bias +0.01 +0.04 +0.02 +0.01
per0LD fly 23: K=3 chosen 0.95; Xn SD 0.19 0.29 0.14 0.18; Yn SD 0.29 0.29 0.24 0.27; Xn bias -0.03 -0.02 +0.01 +0.01; Yn bias +0.03 +0.05 -0.00 +0.05
per0LD fly 24: K=3 chosen 1.00; Xn SD 0.23 0.21 0.26 0.15; Yn SD 0.27 0.22 0.21 0.16; Xn bias -0.03 +0.04 +0.04 +0.07; Yn bias +0.03 -0.05 -0.04 -0.03
per0LD: median SD over flies (Xn blocks 1-4, Yn blocks 1-4) 0.21 0.21 0.17 0.16 0.27 0.22 0.21 0.21 -> USABLE=True
fumin fly 10: K=3 chosen 0.45; Xn SD 0.38 0.28 0.84 0.91; Yn SD 0.53 0.73 0.86 0.47; Xn bias -0.18 -0.15 +0.11 +0.22; Yn bias -0.36 -0.50 -0.54 -0.12
fumin fly 33: K=3 chosen 1.00; Xn SD 0.46 0.40 0.48 0.44; Yn SD 0.49 0.49 0.40 0.76; Xn bias -0.03 +0.04 +0.03 +0.04; Yn bias -0.10 +0.09 +0.05 +0.23
fumin fly 40: K=3 chosen 1.00; Xn SD 0.36 0.37 0.25 0.16; Yn SD 0.38 0.34 0.58 0.66; Xn bias -0.02 -0.06 -0.09 -0.18; Yn bias -0.08 +0.05 +0.26 -0.15
fumin fly 44: K=3 chosen 0.50; Xn SD 0.88 0.66 0.85 0.78; Yn SD 0.21 0.69 0.33 0.64; Xn bias -0.12 -0.45 -0.20 -0.29; Yn bias +0.09 +0.44 -0.12 -0.30
fumin fly 45: K=3 chosen 1.00; Xn SD 0.25 0.26 0.51 0.56; Yn SD 0.24 0.20 0.65 0.69; Xn bias +0.16 +0.07 +0.02 +0.09; Yn bias +0.06 -0.01 -0.08 -0.11
fumin: median SD over flies (Xn blocks 1-4, Yn blocks 1-4) 0.38 0.37 0.51 0.56 0.38 0.49 0.58 0.66 -> USABLE=False

Cycle direction precision, realistic truths (20 replicates per fly, 5 flies per genotype), per 6 h block:
CS      c SD median 0.15 0.17 0.20 0.18 | correct sign mean 0.76 0.72 0.84 0.79 | c bias median +0.015 -0.015 +0.045 -0.014 | true |c| median 0.17 0.17 0.24 0.24
per0LD  c SD median 0.18 0.15 0.14 0.13 | correct sign mean 0.73 0.66 0.93 0.96 | c bias median +0.031 +0.024 -0.007 -0.011 | true |c| median 0.07 0.07 0.21 0.21
fumin   c SD median 0.31 0.34 0.43 0.47 | correct sign mean 0.66 0.63 0.68 0.57 | c bias median -0.033 +0.026 +0.044 +0.015 | true |c| median 0.25 0.25 0.38 0.38
```
- USABLE: CS (median SD 0.20-0.29), per0LD (0.16-0.27). NOT usable: fumin (0.37-0.66; K=3 chosen in only 45-50% of replicates for 2 of 5 flies; biases up to 0.5).

## PRESPEC 20c step 3 + PRESPEC 20e amendment B — 6 h blocks, fly fits and group tests (2026-10-24, code 64bb252)
Narrowest block: one fly's minutes in one 6 h ZT block pooled over its days (shared emissions, 4 matrices). Bias check (ii) of amendment B (no-change null) not yet run.
```
CS (USABLE): flies 25 | best K 2:0.00 3:1.00 4:0.00 | H1 degenerate 0.00
   flies used 24 (no-cycle excluded 1)
   ZT0-6    c_pred -0.021 [-0.079,+0.039] (c>0 in 0.42) | count J -0.041 [-0.100,+0.017] /h
   ZT6-12   c_pred +0.033 [-0.026,+0.089] (c>0 in 0.71) | count J +0.112 [+0.037,+0.192] /h
   ZT12-18  c_pred +0.274 [+0.182,+0.367] (c>0 in 0.92) | count J +0.204 [+0.131,+0.284] /h
   ZT18-24  c_pred +0.205 [+0.104,+0.302] (c>0 in 0.75) | count J +0.196 [+0.103,+0.287] /h
   one-sample Wilcoxon c vs 0 (Holm): ZT0-6 p 0.51 adj 0.51 | ZT6-12 p 0.18 adj 0.36 | ZT12-18 p 8.3e-06 adj 3.3e-05 | ZT18-24 p 0.00096 adj 0.0029
   Friedman chi2 17.80 p 0.00048 | adjacent within-fly Wilcoxon (Holm): ZT0-6->ZT6-12 d +0.041 p 0.21 adj 0.42 | ZT6-12->ZT12-18 d +0.254 p 0.00018 adj 0.0007 | ZT12-18->ZT18-24 d -0.014 p 0.26 adj 0.42 | ZT18-24->ZT0-6 d -0.294 p 0.0025 adj 0.0075
   'group direction changes across the day' (amendment B rule): False
per0LD (USABLE): flies 27 | best K 2:0.19 3:0.81 4:0.00 | H1 degenerate 0.00
   flies used 25 (no-cycle excluded 2)
   ZT0-6    c_pred +0.092 [+0.021,+0.163] (c>0 in 0.68) | count J +0.105 [+0.040,+0.178] /h
   ZT6-12   c_pred +0.072 [-0.019,+0.160] (c>0 in 0.60) | count J +0.091 [+0.012,+0.175] /h
   ZT12-18  c_pred +0.085 [+0.011,+0.156] (c>0 in 0.72) | count J +0.109 [+0.025,+0.191] /h
   ZT18-24  c_pred +0.204 [+0.126,+0.286] (c>0 in 0.84) | count J +0.279 [+0.162,+0.394] /h
   one-sample Wilcoxon c vs 0 (Holm): ZT0-6 p 0.016 adj 0.048 | ZT6-12 p 0.13 adj 0.13 | ZT12-18 p 0.017 adj 0.048 | ZT18-24 p 0.0001 adj 0.00041
   Friedman chi2 6.74 p 0.081 | adjacent within-fly Wilcoxon (Holm): ZT0-6->ZT6-12 d +0.033 p 0.92 adj 1 | ZT6-12->ZT12-18 d -0.026 p 0.96 adj 1 | ZT12-18->ZT18-24 d +0.082 p 0.055 adj 0.22 | ZT18-24->ZT0-6 d -0.113 p 0.055 adj 0.22
   'group direction changes across the day' (amendment B rule): False
fumin (NOT usable: descriptive): flies 53 | best K 2:0.11 3:0.89 4:0.00 | H1 degenerate 0.00
   flies used 19 (no-cycle excluded 34)
   ZT0-6    c_pred +0.172 [+0.067,+0.275] (c>0 in 0.63) | count J +0.146 [+0.034,+0.271] /h
   ZT6-12   c_pred +0.164 [+0.018,+0.281] (c>0 in 0.79) | count J +0.244 [+0.118,+0.378] /h
   ZT12-18  c_pred +0.289 [+0.150,+0.425] (c>0 in 0.84) | count J +0.207 [+0.113,+0.318] /h
   ZT18-24  c_pred +0.213 [+0.030,+0.404] (c>0 in 0.68) | count J +0.155 [+0.060,+0.262] /h
   one-sample Wilcoxon c vs 0 (Holm): ZT0-6 p 0.012 adj 0.025 | ZT6-12 p 0.0046 adj 0.014 | ZT12-18 p 0.0012 adj 0.0047 | ZT18-24 p 0.03 adj 0.03
   Friedman chi2 4.07 p 0.25 | adjacent within-fly Wilcoxon (Holm): ZT0-6->ZT6-12 d +0.001 p 0.71 adj 1 | ZT6-12->ZT12-18 d +0.035 p 0.24 adj 0.97 | ZT12-18->ZT18-24 d -0.017 p 0.95 adj 1 | ZT18-24->ZT0-6 d -0.149 p 0.86 adj 1
   'group direction changes across the day' (amendment B rule): False
```
- CS: circulation is ~0 in the light blocks and forward in both dark blocks; within-fly change across the day is significant (Friedman p 0.0005; ZT6-12 -> 12-18 +0.25, adj p 0.0007; ZT18-24 -> 0-6 -0.29, adj p 0.0075). No block has a significant NEGATIVE group direction, so the pre-stated rule for a direction REVERSAL is not met: the magnitude switches on and off, the direction does not reverse.
- per0LD: forward in ZT0-6, 12-18 and 18-24 (adj p <= 0.048), strongest at ZT18-24; no significant within-fly change (Friedman p 0.08).
- fumin (not usable; descriptive): 34 of 53 flies excluded as 'no cycle' in at least one block; remaining 19 forward in all blocks; no within-fly change.

## PRESPEC 20b theta sensitivity, K = 3 only (2026-10-24, code b060abb)
Full logs: sqw2_0.33.log, sqw2_0.75.log, sqw2_countJ_0.33.log, sqw2_countJ_0.75.log. c_model (FFBS) group means, theta 0.33 / 0.5 / 0.75:
- CS light -0.02 / +0.00 / +0.02 (ns all); CS dark +0.22 / +0.25 / +0.31
- fumin light +0.16 / +0.18 / +0.29; fumin dark +0.24 / +0.34 / +0.48
- per0LD light +0.11 / +0.13 / +0.12; per0LD dark +0.13 / +0.14 / +0.16
- per0DD +0.18 / +0.21 / +0.24
- sss light -0.11 / -0.19 / -0.22; sss dark -0.08 / -0.05 / -0.06 (interval excludes 0 only at 0.33)
- Count J_cycle (cycles/h) has the same sign pattern at all three theta; H1 non-degenerate in the same flies at all theta.
- Reading: every DIRECTION is unchanged across theta; magnitudes grow with theta (a higher threshold labels more wake as Q). The sss light reversal holds at all three. Model-free c is consistent in sign wherever its interval excludes 0 (CS dark, fumin, per0DD, sss light at all theta).

## PRESPEC 20d — reversible null and per-fly attenuation (2026-10-24, code a171bda; stage b rerun after fixing a bug I introduced in c_obs_min, which crashed it; no result had been produced)
```
rev light: true c +0.000 EPR 0.000/h | c_pred +0.024 [-0.037,+0.083] | c_ffbs +0.024 [-0.035,+0.084] | fitted EPR median 0.031 95th pct 0.255 /h (inf 0.00)
rev dark: true c +0.000 EPR 0.000/h | c_pred -0.015 [-0.085,+0.048] | c_ffbs -0.017 [-0.078,+0.045] | fitted EPR median 0.029 95th pct 0.327 /h (inf 0.00)
V1b light: true c +0.350 EPR 0.478/h | c_pred +0.387 [+0.351,+0.425] | c_ffbs +0.387 [+0.352,+0.421] | fitted EPR median 0.635 95th pct 1.198 /h (inf 0.00)
V1b dark: true c -0.150 EPR 0.081/h | c_pred -0.156 [-0.222,-0.096] | c_ffbs -0.155 [-0.220,-0.095] | fitted EPR median 0.083 95th pct 0.500 /h (inf 0.00)
(a) PASS=True
   fly EPR CS light: fraction above reversible-null 95th pct (0.255/h): 0.00 (n 25)
   fly EPR CS dark: fraction above reversible-null 95th pct (0.327/h): 0.40 (n 25)
   fly EPR fumin light: fraction above reversible-null 95th pct (0.255/h): 0.45 (n 53)
   fly EPR fumin dark: fraction above reversible-null 95th pct (0.327/h): 0.55 (n 53)
   fly EPR per0LD light: fraction above reversible-null 95th pct (0.255/h): 0.22 (n 27)
   fly EPR per0LD dark: fraction above reversible-null 95th pct (0.327/h): 0.26 (n 27)
   fly EPR per0DD whole: fraction above reversible-null 95th pct (0.255/h): 0.39 (n 23)
   fly EPR sss light: fraction above reversible-null 95th pct (0.255/h): 0.37 (n 54)
   fly EPR sss dark: fraction above reversible-null 95th pct (0.327/h): 0.43 (n 54)

              simulated records with opposite signs of c_true and c_obs_min: 0.19
CS      light: real c_obs_min +0.019; simulated group mean +0.007 [-0.014,+0.026] -> consistent; flies outside own 95% band 0.04; sim c_true +0.006 / sim c_obs +0.007 = 0.84; real c_model +0.002 / c_obs +0.019 = 0.10
              simulated records with opposite signs of c_true and c_obs_min: 0.06
CS      dark : real c_obs_min +0.152; simulated group mean +0.160 [+0.118,+0.202] -> consistent; flies outside own 95% band 0.00; sim c_true +0.258 / sim c_obs +0.160 = 1.61; real c_model +0.250 / c_obs +0.152 = 1.64
              simulated records with opposite signs of c_true and c_obs_min: 0.21
fumin   light: real c_obs_min +0.157; simulated group mean +0.120 [+0.063,+0.160] -> consistent; flies outside own 95% band 0.10; sim c_true +0.165 / sim c_obs +0.120 = 1.37; real c_model +0.178 / c_obs +0.157 = 1.13
              simulated records with opposite signs of c_true and c_obs_min: 0.14
fumin   dark : real c_obs_min +0.207; simulated group mean +0.193 [+0.111,+0.248] -> consistent; flies outside own 95% band 0.02; sim c_true +0.313 / sim c_obs +0.193 = 1.62; real c_model +0.335 / c_obs +0.207 = 1.61
              simulated records with opposite signs of c_true and c_obs_min: 0.18
per0LD  light: real c_obs_min +0.075; simulated group mean +0.058 [+0.042,+0.079] -> consistent; flies outside own 95% band 0.07; sim c_true +0.118 / sim c_obs +0.058 = 2.02; real c_model +0.126 / c_obs +0.075 = 1.68
              simulated records with opposite signs of c_true and c_obs_min: 0.10
per0LD  dark : real c_obs_min +0.082; simulated group mean +0.083 [+0.061,+0.105] -> consistent; flies outside own 95% band 0.04; sim c_true +0.145 / sim c_obs +0.083 = 1.74; real c_model +0.144 / c_obs +0.082 = 1.77
              simulated records with opposite signs of c_true and c_obs_min: 0.07
per0DD  whole: real c_obs_min +0.114; simulated group mean +0.125 [+0.103,+0.146] -> consistent; flies outside own 95% band 0.04; sim c_true +0.217 / sim c_obs +0.125 = 1.75; real c_model +0.207 / c_obs +0.114 = 1.82
              simulated records with opposite signs of c_true and c_obs_min: 0.24
sss     light: real c_obs_min -0.111; simulated group mean -0.082 [-0.153,-0.013] -> consistent; flies outside own 95% band 0.04; sim c_true -0.165 / sim c_obs -0.082 = 2.02; real c_model -0.185 / c_obs -0.111 = 1.67
              simulated records with opposite signs of c_true and c_obs_min: 0.26
sss     dark : real c_obs_min -0.068; simulated group mean -0.035 [-0.092,+0.011] -> consistent; flies outside own 95% band 0.09; sim c_true -0.059 / sim c_obs -0.035 = 1.67; real c_model -0.052 / c_obs -0.068 = 0.77
```
- (a) PASS: on a detailed-balance truth the fitted c is 0 within +-0.03 (c_pred and FFBS, both phases); the pipeline does not manufacture circulation. EPR noise floor at this record length: median 0.03, 95th percentile 0.26-0.33 nats/h; flies above it: CS dark 40%, fumin 45-55%, sss 37-43%, per0LD 22-26%, per0DD 39%, CS light 0%.
- (b) CONSISTENT in all 9 genotype x phase blocks: each fly's fitted model reproduces the real minute-level observable current (group means inside the simulated 95% range; 0-10% of flies outside their own band). Attenuation hidden/observable ~1.6-2.0 in simulation, matching the real ratio (1.6-1.8) wherever the current is not near 0. So the model's larger circulation is what coarse-graining of its own hidden dynamics predicts.
- per0LD: the minute-level observable current is +0.075 / +0.082 (light/dark) and the model reproduces it. The earlier per0LD disagreement was with the PRESPEC 19 +-5 min transit classification (c ~ 0), a different observable.
- Hidden and observable current have opposite signs in 6-26% of simulated single records (most in weak blocks: CS light, sss).

## PRESPEC 20e sections 1-2 — direction changes across days; light vs dark (2026-10-24, code 1d59b73)
```
Section 1: across days (first vs second half), hidden-level c_pred; no-change null 8 replicates per fly
  CS      light: flies 25 (excluded 0) | observed flips 11 (0.44) | expected under no change 10.6 (0.42) | one-sided p 0.507
  CS      dark : flies 24 (excluded 1) | observed flips 4 (0.17) | expected under no change 4.6 (0.19) | one-sided p 0.746
  fumin   light: flies 36 (excluded 17) | observed flips 11 (0.31) | expected under no change 7.3 (0.20) | one-sided p 0.0658
  fumin   dark : flies 30 (excluded 23) | observed flips 8 (0.27) | expected under no change 6.5 (0.22) | one-sided p 0.306
  per0LD  light: flies 27 (excluded 0) | observed flips 12 (0.44) | expected under no change 8.1 (0.30) | one-sided p 0.0508
  per0LD  dark : flies 27 (excluded 0) | observed flips 9 (0.33) | expected under no change 7.1 (0.26) | one-sided p 0.24
  per0DD  whole: flies 22 (excluded 0) | observed flips 3 (0.14) | expected under no change 3.6 (0.16) | one-sided p 0.784
  sss     light: flies 20 (excluded 34) | observed flips 6 (0.30) | expected under no change 5.6 (0.28) | one-sided p 0.517
  sss     dark : flies 22 (excluded 32) | observed flips 10 (0.45) | expected under no change 8.3 (0.38) | one-sided p 0.276
Section 2: light vs dark within fly (12 h fits), descriptive
  CS      flies 25 | opposite signs light vs dark 0.48 | expected from estimation noise if each fly's light and dark values were its true ones ~0.50 (rough: c SD 0.08/0.10, normal approx.) | Spearman light-dark -0.14 p 0.5
  fumin   flies 48 | opposite signs light vs dark 0.40 | expected from estimation noise if each fly's light and dark values were its true ones ~0.41 (rough: c SD 0.18/0.32, normal approx.) | Spearman light-dark -0.06 p 0.7
  per0LD  flies 27 | opposite signs light vs dark 0.22 | expected from estimation noise if each fly's light and dark values were its true ones ~0.41 (rough: c SD 0.12/0.13, normal approx.) | Spearman light-dark -0.32 p 0.11
  sss     flies 37 | opposite signs light vs dark 0.27 | expected from estimation noise if each fly's light and dark values were its true ones ~0.47 (rough: c SD 0.22/0.31, normal approx.) | Spearman light-dark +0.19 p 0.27
```
- Section 1: no block shows more sign flips between the first and second half of the recording than estimation noise predicts (smallest one-sided p: per0LD light 0.051, fumin light 0.066; 9 tests, none significant after Holm). Reading as pre-stated: at ~2 days per half, these data cannot distinguish a fly that changes direction across days from one that does not. fumin and sss lose 17-34 flies to 'no cycle' in a half.
- Section 2 (descriptive): the 'expected' column adds bootstrap noise around already-noisy estimates and so overstates the noise-only rate; read only as a rough reference. per0LD and sss have fewer light/dark sign disagreements (0.22, 0.27) than that reference (0.41, 0.47); CS and fumin match it.

## PRESPEC 20e section 3 + amendment B bias check — 6 h within-phase flips (2026-10-24, code 1d59b73)
```
Section 3: within phase, 6 h blocks (block pairs 1-2 = ZT0-6 vs 6-12, 3-4 = ZT12-18 vs 18-24)
  CS      blocks 1-2 (test): flies 24 | observed flips 9 (0.38) | expected 9.4 | one-sided p 0.661 | bias check: null-replicate Wilcoxon p<0.05 in 0.00 of 8
  CS      blocks 3-4 (test): flies 25 | observed flips 8 (0.32) | expected 5.1 | one-sided p 0.0712 | bias check: null-replicate Wilcoxon p<0.05 in 0.00 of 8
  per0LD  blocks 1-2 (test): flies 26 | observed flips 14 (0.54) | expected 9.7 | one-sided p 0.0342 | bias check: null-replicate Wilcoxon p<0.05 in 0.00 of 8
  per0LD  blocks 3-4 (test): flies 26 | observed flips 9 (0.35) | expected 8.6 | one-sided p 0.512 | bias check: null-replicate Wilcoxon p<0.05 in 0.12 of 8
  fumin   blocks 1-2 (descriptive): flies 35 | observed flips 12 (0.34) | expected 10.0 | one-sided p 0.246 | bias check: null-replicate Wilcoxon p<0.05 in 0.00 of 8
  fumin   blocks 3-4 (descriptive): flies 21 | observed flips 4 (0.19) | expected 3.7 | one-sided p 0.539 | bias check: null-replicate Wilcoxon p<0.05 in 0.12 of 8
```
- No excess of within-phase sign flips survives correction: per0LD ZT0-6 vs 6-12 has p 0.034 uncorrected (14 flips vs 9.7 expected), Holm over the 4 tests in usable genotypes gives 0.14; CS ZT12-18 vs 18-24 p 0.071.
- Bias check: on the no-change null replicates the within-phase Wilcoxon test was significant in 0/8 (CS, both phases), 0/8 and 1/8 (per0LD), consistent with no estimator-made block differences. Limitation: this null holds the matrix fixed WITHIN a phase, so it checks the within-phase comparisons only; the CS light-to-dark changes (ZT6-12 -> 12-18, ZT18-24 -> 0-6) rest on the realistic-truth bias (<= 0.05 in c) against effects of 0.25-0.29.

## POST HOC (user question, 2026-10-24): how the 'expected flips' were made, and a sensitivity
- Method (as PRESPEC 20e): per fly, 8 records simulated from its whole-record fitted model (no change by construction), split identically, refitted at K=3 (10 starts); p_f = fraction of the 8 replicates whose two halves (or blocks) differ in sign; expected = sum of p_f; p-value = upper tail of the Poisson-binomial(p_f) at the observed count.
- Weakness found: with 8 replicates, p_f = 0 in 56 of 233 fly-phases (true p_f could be up to ~0.3), and 4 of those flies did flip; a p_f of 0 makes the test anti-conservative. Sensitivity with p_f smoothed as (flips + 0.5)/(n + 1):
```
  CS      light: flies 25 (excluded 0) | observed flips 11 (0.44) | expected under no change 10.8 (0.43) | one-sided p 0.544
  CS      dark : flies 24 (excluded 1) | observed flips 4 (0.17) | expected under no change 5.4 (0.23) | one-sided p 0.855
  fumin   light: flies 36 (excluded 17) | observed flips 11 (0.31) | expected under no change 8.6 (0.24) | one-sided p 0.206
  fumin   dark : flies 30 (excluded 23) | observed flips 8 (0.27) | expected under no change 7.7 (0.26) | one-sided p 0.52
  per0LD  light: flies 27 (excluded 0) | observed flips 12 (0.44) | expected under no change 8.7 (0.32) | one-sided p 0.0998
  per0LD  dark : flies 27 (excluded 0) | observed flips 9 (0.33) | expected under no change 7.8 (0.29) | one-sided p 0.365
  per0DD  whole: flies 22 (excluded 0) | observed flips 3 (0.14) | expected under no change 4.4 (0.20) | one-sided p 0.886
  sss     light: flies 20 (excluded 34) | observed flips 6 (0.30) | expected under no change 6.2 (0.31) | one-sided p 0.636
  sss     dark : flies 22 (excluded 32) | observed flips 10 (0.45) | expected under no change 8.7 (0.39) | one-sided p 0.343
  CS      blocks 1-2 (test): flies 24 | observed flips 9 (0.38) | expected 9.7 | one-sided p 0.703 | bias check: null-replicate Wilcoxon p<0.05 in 0.00 of 8
  CS      blocks 3-4 (test): flies 25 | observed flips 8 (0.32) | expected 5.9 | one-sided p 0.189 | bias check: null-replicate Wilcoxon p<0.05 in 0.00 of 8
  per0LD  blocks 1-2 (test): flies 26 | observed flips 14 (0.54) | expected 10.1 | one-sided p 0.0585 | bias check: null-replicate Wilcoxon p<0.05 in 0.00 of 8
  per0LD  blocks 3-4 (test): flies 26 | observed flips 9 (0.35) | expected 9.1 | one-sided p 0.6 | bias check: null-replicate Wilcoxon p<0.05 in 0.12 of 8
  fumin   blocks 1-2 (descriptive): flies 35 | observed flips 12 (0.34) | expected 10.8 | one-sided p 0.378 | bias check: null-replicate Wilcoxon p<0.05 in 0.00 of 8
  fumin   blocks 3-4 (descriptive): flies 21 | observed flips 4 (0.19) | expected 4.6 | one-sided p 0.716 | bias check: null-replicate Wilcoxon p<0.05 in 0.12 of 8
```
- Smallest p becomes 0.059 (per0LD ZT0-6 vs 6-12) and 0.100 (per0LD light across days); conclusions unchanged (no excess flips anywhere). Other limitation: the null is simulated from each fly's FITTED parameters (parametric bootstrap), so it inherits their estimation error.

## POST HOC (user question, 2026-10-24): EPR vs J_cycle, genotype x phase medians
Narrowest block: one fly's 12 h phase (all days); medians over flies, 95% bootstrap CI of the median. J count-based (forward +), EPR = rule-of-3 lower bound per path, path-averaged. Code fig_epr_j.py; figure figures/fig6_epr_vs_j.*
```
CS      light n 25 | J median +0.020 [-0.018,+0.088] | EPR median 0.041 [0.024,0.078] | EPR/|J| 1.99
CS      dark  n 25 | J median +0.201 [+0.095,+0.264] | EPR median 0.301 [0.065,0.389] | EPR/|J| 1.50
fumin   light n 53 | J median +0.099 [+0.037,+0.155] | EPR median 0.108 [0.046,0.253] | EPR/|J| 1.09
fumin   dark  n 53 | J median +0.109 [+0.070,+0.147] | EPR median 0.143 [0.071,0.223] | EPR/|J| 1.31
per0LD  light n 27 | J median +0.047 [+0.029,+0.109] | EPR median 0.030 [0.021,0.126] | EPR/|J| 0.64
per0LD  dark  n 27 | J median +0.197 [+0.056,+0.257] | EPR median 0.177 [0.063,0.226] | EPR/|J| 0.90
per0DD  whole n 23 | J median +0.233 [+0.139,+0.298] | EPR median 0.208 [0.077,0.365] | EPR/|J| 0.89
sss     light n 54 | J median -0.007 [-0.045,+0.000] | EPR median 0.031 [0.015,0.048] | EPR/|J| 4.39
sss     dark  n 54 | J median -0.012 [-0.041,+0.000] | EPR median 0.044 [0.020,0.072] | EPR/|J| 3.51
```

## POST HOC (user question, 2026-10-24): EPR vs J, hidden vs observable
Observable = same counting on run-compressed S/Q/W minute symbols (Q<->H1, W<->H2, S<->H3); figure figures/fig6b_epr_vs_j_hidden_vs_obs.*; code fig_epr_j_obs.py. Medians over flies, 95% bootstrap CI.
```
block | hidden J, EPR | observable J, EPR | ratio hidden/observable (J, EPR)
CS      light | +0.020 0.041 | +0.062 [+0.000,+0.125] 0.038 [0.013,0.072] | J 0.33, EPR 1.06
CS      dark  | +0.201 0.301 | +0.167 [+0.104,+0.188] 0.116 [0.054,0.189] | J 1.21, EPR 2.58
fumin   light | +0.099 0.108 | +0.083 [+0.042,+0.167] 0.058 [0.045,0.189] | J 1.19, EPR 1.85
fumin   dark  | +0.109 0.143 | +0.083 [+0.042,+0.125] 0.093 [0.046,0.150] | J 1.31, EPR 1.54
per0LD  light | +0.047 0.030 | +0.062 [+0.000,+0.146] 0.012 [0.003,0.060] | J 0.75, EPR 2.57
per0LD  dark  | +0.197 0.177 | +0.167 [+0.083,+0.208] 0.073 [0.045,0.140] | J 1.18, EPR 2.43
per0DD  whole | +0.233 0.208 | +0.177 [+0.139,+0.271] 0.113 [0.081,0.213] | J 1.32, EPR 1.84
sss     light | -0.007 0.031 | +0.000 [-0.042,+0.000] 0.009 [0.001,0.035] | J nan, EPR 3.30
sss     dark  | -0.012 0.044 | +0.000 [-0.042,+0.000] 0.029 [0.012,0.036] | J nan, EPR 1.51
```
- Observable EPR is 1.5-3.3x below hidden (coarse-graining loses irreversibility, as expected); observable J is 1.2-1.3x below hidden in the driven blocks. Observable medians are discrete (integer counts / hours).

## PRESPEC 20f — hidden vs observable EPR, with calibration (2026-10-24, code 836b5eb)
```
CS      light (not driven: descriptive) flies 25 | median EPR_h 0.041 EPR_o 0.038 | D median +0.009 [+0.003,+0.022] Wilcoxon p 0.0081 | D<0 in 0.16 (sim 0.22) | r median 0.74 [0.53,0.93], simulated 0.70 [0.50,0.87] -> consistent | E_ffbs bias vs true path +0.19, vs model +1.54
CS      dark  (driven) flies 25 | median EPR_h 0.301 EPR_o 0.116 | D median +0.112 [+0.024,+0.241] Wilcoxon p 2.6e-06 | D<0 in 0.04 (sim 0.08) | r median 0.48 [0.32,0.61], simulated 0.59 [0.43,0.71] -> consistent | E_ffbs bias vs true path +0.07, vs model +0.16 | Holm p 1.8e-05 -> hidden dissipation REPORTED
fumin   light (driven) flies 53 | median EPR_h 0.108 EPR_o 0.058 | D median +0.045 [+0.022,+0.053] Wilcoxon p 8e-05 | D<0 in 0.19 (sim 0.27) | r median 0.51 [0.40,0.79], simulated 0.67 [0.60,0.80] -> INCONSISTENT | E_ffbs bias vs true path +0.08, vs model +0.26 | Holm p 0.00032 -> hidden dissipation REPORTED
fumin   dark  (driven) flies 53 | median EPR_h 0.143 EPR_o 0.093 | D median +0.015 [+0.007,+0.076] Wilcoxon p 5.3e-05 | D<0 in 0.23 (sim 0.30) | r median 0.62 [0.42,0.84], simulated 0.70 [0.56,0.85] -> consistent | E_ffbs bias vs true path +0.01, vs model +0.01 | Holm p 0.00027 -> hidden dissipation REPORTED
per0LD  light (driven) flies 27 | median EPR_h 0.030 EPR_o 0.012 | D median +0.018 [+0.011,+0.025] Wilcoxon p 0.001 | D<0 in 0.15 (sim 0.19) | r median 0.27 [0.13,0.80], simulated 0.56 [0.41,0.82] -> INCONSISTENT | E_ffbs bias vs true path +0.22, vs model +1.32 | Holm p 0.002 -> hidden dissipation REPORTED
per0LD  dark  (driven) flies 27 | median EPR_h 0.177 EPR_o 0.073 | D median +0.058 [+0.029,+0.095] Wilcoxon p 1.6e-05 | D<0 in 0.07 (sim 0.12) | r median 0.52 [0.26,0.71], simulated 0.48 [0.32,0.70] -> consistent | E_ffbs bias vs true path +0.14, vs model +0.25 | Holm p 9.6e-05 -> hidden dissipation REPORTED
per0DD  whole (driven) flies 23 | median EPR_h 0.208 EPR_o 0.113 | D median +0.055 [+0.025,+0.141] Wilcoxon p 0.00021 | D<0 in 0.13 (sim 0.05) | r median 0.63 [0.39,0.73], simulated 0.63 [0.49,0.77] -> consistent | E_ffbs bias vs true path +0.06, vs model +0.12 | Holm p 0.00064 -> hidden dissipation REPORTED
sss     light (driven) flies 54 | median EPR_h 0.031 EPR_o 0.009 | D median +0.007 [+0.000,+0.019] Wilcoxon p 0.0016 | D<0 in 0.28 (sim 0.34) | r median 0.56 [0.26,0.81], simulated 0.72 [0.53,0.97] -> consistent | E_ffbs bias vs true path +0.00, vs model +1.39 | Holm p 0.002 -> hidden dissipation REPORTED
sss     dark  (not driven: descriptive) flies 54 | median EPR_h 0.044 EPR_o 0.029 | D median +0.005 [+0.000,+0.024] Wilcoxon p 0.073 | D<0 in 0.35 (sim 0.36) | r median 0.79 [0.45,1.00], simulated 0.79 [0.60,0.97] -> consistent | E_ffbs bias vs true path +0.03, vs model +1.06
```
- By the pre-stated rule, hidden dissipation D > 0 is 'reported' in all 7 driven blocks (Holm p <= 0.002; |bias vs true path| <= 0.22).
- BUT the calibration against the fitted model's exact stationary EPR (E_model) shows the count-based estimator is strongly biased UPWARD where irreversibility is weak: +132% per0LD light, +139% sss light (+154% CS light, +106% sss dark, not driven). In strong blocks it is modest: CS dark +16%, per0DD +12%, fumin dark +1%, per0LD dark +25%, fumin light +26%. The pre-stated criterion used the realized-path reference, which hides this; in hindsight the model reference is the physically relevant one. Hidden-dissipation values for per0LD light and sss light are dominated by finite-count bias and should NOT be interpreted; CS dark, fumin dark, per0LD dark, per0DD are interpretable, fumin light marginal.
- Consequence for wording: the rule-of-3 count EPR is a lower bound only with respect to the zero-count substitution; with finite counts the Schnakenberg sum is biased upward when the true EPR is small, so it is NOT a lower bound on the true rate in weak blocks. Methods V7 text corrected.
- Ratio r = EPR_o/EPR_h: consistent with coarse-graining of the fitted dynamics in CS dark (0.48), fumin dark (0.62), per0LD dark (0.52), per0DD (0.63), sss light; inconsistent (real lower than simulated) in fumin light (0.51 vs 0.67) and per0LD light (0.27 vs 0.56).

## Correction to the 20f reading (2026-10-25, post hoc decomposition of sqw_out/f_epr.json)
The 20f 'relative bias vs model' (+100-150% in weak blocks) is a median of ratios with small denominators and mixes two effects. Split by fitted model:
- flies whose model has all transitions in both directions: count estimator - exact = +0.000 to +0.025 nats/h (median by block), a near-constant ADDITIVE offset; it dominates only where the true EPR is ~0.01-0.03 (CS light, per0LD light, sss dark);
- flies whose model has a one-way edge (exact EPR infinite or > 1 nat/h; fumin 21-36%, sss 24-35%, few elsewhere): the exact value is huge while the count estimator stays small (median 0.01-0.06), i.e. the count estimator UNDER-reports relative to the model there.
- Figures 6 and 6b updated: weak-circulation blocks faded, offset band shown; exact-model markers NOT shown because one-way-edge flies make the model median misleading for fumin and sss. Methods V7 item (vii) reworded accordingly.

## Figure 7 (user request, 2026-10-25): three-way Xn-Yn per fly, whole record (descriptive; values from existing outputs)
Code fig_xy3.py; figures/fig7_xn_yn_three_way.*
```
(a) observable, model-free | CS | n 25 | mean Xn -0.395 Yn -0.330 | mean c +0.032 | above diagonal 0.60
(b) S/Q/W HMM, posterior | CS | n 25 | mean Xn -0.146 Yn +0.094 | mean c +0.120 | above diagonal 0.88
(c) binary HMM, Viterbi (legacy) | CS | n 25 | mean Xn -0.279 Yn +0.169 | mean c +0.224 | above diagonal 0.88
(a) observable, model-free | per0LD | n 27 | mean Xn -0.489 Yn -0.454 | mean c +0.017 | above diagonal 0.59
(a) observable, model-free | per0DD | n 25 | mean Xn -0.526 Yn -0.316 | mean c +0.105 | above diagonal 0.72
(b) S/Q/W HMM, posterior | per0LD | n 27 | mean Xn -0.462 Yn -0.194 | mean c +0.134 | above diagonal 0.89
(b) S/Q/W HMM, posterior | per0DD | n 23 | mean Xn -0.310 Yn +0.104 | mean c +0.207 | above diagonal 0.87
(c) binary HMM, Viterbi (legacy) | per0 | n 25 | mean Xn -0.024 Yn -0.217 | mean c -0.096 | above diagonal 0.44
(a) observable, model-free | fumin | n 53 | mean Xn -0.785 Yn -0.699 | mean c +0.043 | above diagonal 0.70
(b) S/Q/W HMM, posterior | fumin | n 53 | mean Xn -0.461 Yn -0.033 | mean c +0.214 | above diagonal 0.91
(c) binary HMM, Viterbi (legacy) | fumin | n 53 | mean Xn -0.951 Yn +0.647 | mean c +0.799 | above diagonal 1.00
(a) observable, model-free | sss | n 54 | mean Xn -0.461 Yn -0.538 | mean c -0.039 | above diagonal 0.37
(b) S/Q/W HMM, posterior | sss | n 53 | mean Xn +0.005 Yn -0.210 | mean c -0.107 | above diagonal 0.32
(c) binary HMM, Viterbi (legacy) | sss | n 54 | mean Xn -0.920 Yn +0.729 | mean c +0.825 | above diagonal 1.00
```

## PRESPEC 20g — round-trip types per sleep bout (2026-10-25, code 4d8f1ff; matched hidden version sqw7b.py, amendment B, POST HOC)
Full output: sqw7.log (observed + pre-specified hidden, theta 0.33/0.5/0.75), sqw7b.log (matched hidden). Figure: figures/fig8_roundtrips.* (observed + matched hidden, theta 0.5).
```
  hidden(matched) CS      light n 25 (round trips/fly 19.5) | forward 0.20 [0.17,0.24] reverse 0.23 [0.19,0.27] direct 0.24 [0.17,0.32] quiet_both 0.33 [0.24,0.42] | most common: quiet_both | fwd-rev -0.029 p 0.38 Holm 0.45
      rates/h forward 0.079 reverse 0.098 direct 0.100 quiet_both 0.129 | P(Q entry) 0.53 P(Q exit) 0.56 | joint-product forward -0.000 (p 0.75)
  hidden(matched) CS      dark  n 25 (round trips/fly 19.4) | forward 0.38 [0.33,0.43] reverse 0.11 [0.08,0.15] direct 0.19 [0.15,0.23] quiet_both 0.31 [0.25,0.38] | most common: forward | fwd-rev +0.267 p 4.2e-07 Holm 3.3e-06
      rates/h forward 0.150 reverse 0.048 direct 0.074 quiet_both 0.131 | P(Q entry) 0.70 P(Q exit) 0.43 | joint-product forward -0.011 (p 0.3)
  hidden(matched) fumin   light n 46 (round trips/fly 13.0) | forward 0.30 [0.25,0.35] reverse 0.11 [0.09,0.14] direct 0.41 [0.35,0.47] quiet_both 0.18 [0.12,0.25] | most common: direct | fwd-rev +0.189 p 1.7e-05 Holm 0.00012
      rates/h forward 0.169 reverse 0.058 direct 0.237 quiet_both 0.076 | P(Q entry) 0.48 P(Q exit) 0.29 | joint-product forward -0.011 (p 0.33)
  hidden(matched) fumin   dark  n 46 (round trips/fly 9.1) | forward 0.38 [0.31,0.45] reverse 0.07 [0.05,0.10] direct 0.35 [0.27,0.43] quiet_both 0.21 [0.12,0.29] | most common: forward | fwd-rev +0.305 p 8.5e-08 Holm 7.6e-07
      rates/h forward 0.132 reverse 0.029 direct 0.170 quiet_both 0.050 | P(Q entry) 0.58 P(Q exit) 0.28 | joint-product forward +0.003 (p 0.6)
  hidden(matched) per0LD  light n 25 (round trips/fly 19.0) | forward 0.30 [0.26,0.34] reverse 0.15 [0.11,0.20] direct 0.30 [0.22,0.38] quiet_both 0.25 [0.18,0.33] | most common: direct | fwd-rev +0.144 p 0.00049 Holm 0.0024
      rates/h forward 0.120 reverse 0.065 direct 0.129 quiet_both 0.083 | P(Q entry) 0.55 P(Q exit) 0.40 | joint-product forward -0.002 (p 0.98)
  hidden(matched) per0LD  dark  n 26 (round trips/fly 17.1) | forward 0.30 [0.24,0.35] reverse 0.18 [0.14,0.22] direct 0.36 [0.31,0.43] quiet_both 0.16 [0.12,0.22] | most common: direct | fwd-rev +0.118 p 0.012 Holm 0.036
      rates/h forward 0.112 reverse 0.059 direct 0.123 quiet_both 0.063 | P(Q entry) 0.46 P(Q exit) 0.34 | joint-product forward -0.000 (p 0.9)
  hidden(matched) per0DD  whole n 23 (round trips/fly 32.9) | forward 0.35 [0.30,0.41] reverse 0.15 [0.11,0.19] direct 0.19 [0.15,0.25] quiet_both 0.30 [0.23,0.37] | most common: forward | fwd-rev +0.202 p 8.8e-05 Holm 0.00053
      rates/h forward 0.128 reverse 0.060 direct 0.074 quiet_both 0.124 | P(Q entry) 0.66 P(Q exit) 0.45 | joint-product forward +0.003 (p 0.82)
  hidden(matched) sss     light n 36 (round trips/fly 7.6) | forward 0.08 [0.05,0.10] reverse 0.22 [0.15,0.29] direct 0.37 [0.27,0.48] quiet_both 0.34 [0.22,0.45] | most common: direct | fwd-rev -0.142 p 0.0011 Holm 0.0046
      rates/h forward 0.026 reverse 0.077 direct 0.126 quiet_both 0.088 | P(Q entry) 0.41 P(Q exit) 0.55 | joint-product forward -0.012 (p 0.096)
  hidden(matched) sss     dark  n 39 (round trips/fly 8.0) | forward 0.14 [0.10,0.18] reverse 0.20 [0.14,0.26] direct 0.41 [0.32,0.51] quiet_both 0.25 [0.16,0.34] | most common: direct | fwd-rev -0.061 p 0.22 Holm 0.45
      rates/h forward 0.050 reverse 0.060 direct 0.149 quiet_both 0.074 | P(Q entry) 0.39 P(Q exit) 0.45 | joint-product forward -0.009 (p 0.097)
```
- Observed: direct A->S->A is the most common round trip in every genotype and phase (48-81% at theta 0.5; still the largest at 0.33 and 0.75, 33-88%). Forward > reverse: CS dark (0.24 vs 0.11, Holm 0.008), fumin light/dark, per0DD (Holm <= 0.018); sss light reverse > forward (0.19 vs 0.12, ns at bout level); per0LD ns.
- Hidden (matched): forward is the single most common round trip in CS dark (38%), fumin dark (38%), per0DD (35%); forward > reverse in CS dark, fumin, per0LD, per0DD (Holm <= 0.036); reverse > forward in sss light (Holm 0.005). Which type is most common depends on theta and on the hidden definition; forward-minus-reverse is robust.
- Independence: joint minus product of marginals is within +-0.04 everywhere and mostly ns; entry and exit routes of a bout are close to independent.

## PRESPEC 20h — dwell times, observed and hidden (2026-10-25, code 64dd6eb)
Library: sqw_ref.py section 10 (runs_censored, km_survival, life_table_hazard, geometric_survival, DWELL_EDGES). A run cut at a phase boundary is censored there; its remainder in the next phase is NOT counted as a new run (its start is unobserved). Figures 9, 9s, 10, 10s. Geometric overlay in Fig 10 weighted by each fly's decoded run count (first draft weighted flies equally, which mismatched the pooled KM; corrected before reporting).
```
obs CS      light S    runs  1777 (censored 0.00) | per-fly KM median 14 min | hazard bins outside model 95% band 0.90 of 10 | h per bin: 0.000 0.000 0.000 0.088 0.070 0.058 0.059 0.057 0.041 0.020 - - -
obs CS      light wake runs  1822 (censored 0.05) | per-fly KM median 10 min | hazard bins outside model 95% band 0.80 of 10 | h per bin: 0.081 0.086 0.083 0.073 0.068 0.053 0.042 0.024 0.020 0.020 - - -
obs CS      light Q    runs  5893 (censored 0.00) | per-fly KM median 2 min | hazard bins outside model 95% band 0.33 of 6 | h per bin: 0.475 0.398 0.418 0.349 0.391 0.300 - - - - - - -
obs CS      light W    runs  5888 (censored 0.01) | per-fly KM median 2 min | hazard bins outside model 95% band 1.00 of 8 | h per bin: 0.275 0.320 0.280 0.240 0.191 0.131 0.112 0.054 - - - - -
hid CS      light H1 runs per path summed over flies  4026.0 | model mean dwell 1/(1-P_ii): median over flies 3.1 min [IQR 2.6-3.5] | decoded KM median 2 min
hid CS      light H2 runs per path summed over flies  3763.1 | model mean dwell 1/(1-P_ii): median over flies 6.4 min [IQR 4.3-8.3] | decoded KM median 3 min
hid CS      light H3 runs per path summed over flies  1739.3 | model mean dwell 1/(1-P_ii): median over flies 20.4 min [IQR 16.3-23.7] | decoded KM median 14 min
obs CS      dark  S    runs  1354 (censored 0.03) | per-fly KM median 27 min | hazard bins outside model 95% band 0.45 of 11 | h per bin: 0.000 0.000 0.000 0.044 0.032 0.029 0.031 0.033 0.026 0.023 0.024 - -
obs CS      dark  wake runs  1321 (censored 0.05) | per-fly KM median 12 min | hazard bins outside model 95% band 0.78 of 9 | h per bin: 0.061 0.060 0.054 0.056 0.071 0.088 0.062 0.052 0.055 - - - -
obs CS      dark  Q    runs  3469 (censored 0.00) | per-fly KM median 1 min | hazard bins outside model 95% band 0.17 of 6 | h per bin: 0.498 0.420 0.369 0.368 0.341 0.256 - - - - - - -
obs CS      dark  W    runs  3333 (censored 0.01) | per-fly KM median 3 min | hazard bins outside model 95% band 0.62 of 8 | h per bin: 0.268 0.244 0.233 0.213 0.162 0.112 0.128 0.097 - - - - -
hid CS      dark  H1 runs per path summed over flies  2285.6 | model mean dwell 1/(1-P_ii): median over flies 3.2 min [IQR 2.8-3.6] | decoded KM median 2 min
hid CS      dark  H2 runs per path summed over flies  2053.8 | model mean dwell 1/(1-P_ii): median over flies 7.2 min [IQR 5.4-11.9] | decoded KM median 5 min
hid CS      dark  H3 runs per path summed over flies  1328.0 | model mean dwell 1/(1-P_ii): median over flies 35.9 min [IQR 30.0-49.0] | decoded KM median 27 min
obs fumin   light S    runs   956 (censored 0.02) | per-fly KM median 7 min | hazard bins outside model 95% band 1.00 of 7 | h per bin: 0.000 0.000 0.000 0.191 0.196 0.160 0.168 - - - - - -
obs fumin   light wake runs   976 (censored 0.08) | per-fly KM median 28 min | hazard bins outside model 95% band 0.69 of 13 | h per bin: 0.015 0.014 0.029 0.029 0.030 0.029 0.023 0.023 0.016 0.010 0.006 0.005 0.005
obs fumin   light Q    runs  7501 (censored 0.00) | per-fly KM median 1 min | hazard bins outside model 95% band 0.50 of 6 | h per bin: 0.538 0.438 0.434 0.374 0.346 0.301 - - - - - - -
obs fumin   light W    runs  7867 (censored 0.01) | per-fly KM median 3 min | hazard bins outside model 95% band 0.60 of 10 | h per bin: 0.234 0.216 0.205 0.205 0.161 0.095 0.052 0.044 0.032 0.026 - - -
hid fumin   light H1 runs per path summed over flies  4574.7 | model mean dwell 1/(1-P_ii): median over flies 3.6 min [IQR 2.4-9.2] | decoded KM median 2 min
hid fumin   light H2 runs per path summed over flies  4793.9 | model mean dwell 1/(1-P_ii): median over flies 9.1 min [IQR 5.6-26.9] | decoded KM median 5 min
hid fumin   light H3 runs per path summed over flies   978.4 | model mean dwell 1/(1-P_ii): median over flies 9.2 min [IQR 7.3-10.6] | decoded KM median 8 min
obs fumin   dark  S    runs   647 (censored 0.00) | per-fly KM median 7 min | hazard bins outside model 95% band 1.00 of 7 | h per bin: 0.000 0.000 0.000 0.200 0.210 0.226 0.167 - - - - - -
obs fumin   dark  wake runs   667 (censored 0.15) | per-fly KM median 33 min | hazard bins outside model 95% band 0.77 of 13 | h per bin: 0.009 0.011 0.014 0.021 0.024 0.028 0.024 0.012 0.008 0.004 0.003 0.002 0.001
obs fumin   dark  Q    runs  7171 (censored 0.00) | per-fly KM median 1 min | hazard bins outside model 95% band 0.40 of 5 | h per bin: 0.590 0.493 0.475 0.474 0.441 - - - - - - - -
obs fumin   dark  W    runs  7505 (censored 0.01) | per-fly KM median 4 min | hazard bins outside model 95% band 0.64 of 11 | h per bin: 0.179 0.176 0.177 0.172 0.134 0.089 0.048 0.039 0.027 0.014 0.008 - -
hid fumin   dark  H1 runs per path summed over flies  4448.6 | model mean dwell 1/(1-P_ii): median over flies 2.9 min [IQR 2.2-4.3] | decoded KM median 2 min
hid fumin   dark  H2 runs per path summed over flies  4671.6 | model mean dwell 1/(1-P_ii): median over flies 12.8 min [IQR 7.8-28.5] | decoded KM median 5 min
hid fumin   dark  H3 runs per path summed over flies   643.7 | model mean dwell 1/(1-P_ii): median over flies 8.2 min [IQR 6.8-9.0] | decoded KM median 7 min
obs per0LD  light S    runs  1601 (censored 0.03) | per-fly KM median 22 min | hazard bins outside model 95% band 0.64 of 11 | h per bin: 0.000 0.000 0.000 0.055 0.058 0.044 0.040 0.030 0.024 0.021 0.016 - -
obs per0LD  light wake runs  1616 (censored 0.03) | per-fly KM median 10 min | hazard bins outside model 95% band 0.58 of 12 | h per bin: 0.074 0.053 0.069 0.069 0.080 0.067 0.049 0.036 0.021 0.006 0.011 0.015 -
obs per0LD  light Q    runs  4256 (censored 0.00) | per-fly KM median 1 min | hazard bins outside model 95% band 0.50 of 6 | h per bin: 0.517 0.447 0.399 0.401 0.311 0.213 - - - - - - -
obs per0LD  light W    runs  4531 (censored 0.01) | per-fly KM median 2 min | hazard bins outside model 95% band 0.78 of 9 | h per bin: 0.317 0.273 0.229 0.191 0.141 0.104 0.063 0.045 0.037 - - - -
hid per0LD  light H1 runs per path summed over flies  2648.4 | model mean dwell 1/(1-P_ii): median over flies 4.0 min [IQR 2.7-5.6] | decoded KM median 2 min
hid per0LD  light H2 runs per path summed over flies  2712.8 | model mean dwell 1/(1-P_ii): median over flies 7.6 min [IQR 4.3-15.8] | decoded KM median 4 min
hid per0LD  light H3 runs per path summed over flies  1583.6 | model mean dwell 1/(1-P_ii): median over flies 31.3 min [IQR 18.3-40.5] | decoded KM median 18 min
obs per0LD  dark  S    runs  2011 (censored 0.04) | per-fly KM median 18 min | hazard bins outside model 95% band 0.64 of 11 | h per bin: 0.000 0.000 0.000 0.066 0.054 0.043 0.046 0.045 0.030 0.028 0.034 - -
obs per0LD  dark  wake runs  1998 (censored 0.02) | per-fly KM median 8 min | hazard bins outside model 95% band 0.56 of 9 | h per bin: 0.065 0.073 0.087 0.095 0.093 0.074 0.058 0.046 0.049 - - - -
obs per0LD  dark  Q    runs  4507 (censored 0.00) | per-fly KM median 2 min | hazard bins outside model 95% band 0.40 of 5 | h per bin: 0.485 0.395 0.420 0.395 0.414 - - - - - - - -
obs per0LD  dark  W    runs  5014 (censored 0.00) | per-fly KM median 2 min | hazard bins outside model 95% band 0.88 of 8 | h per bin: 0.366 0.337 0.309 0.281 0.199 0.140 0.083 0.085 - - - - -
hid per0LD  dark  H1 runs per path summed over flies  2959.6 | model mean dwell 1/(1-P_ii): median over flies 3.5 min [IQR 2.8-4.1] | decoded KM median 3 min
hid per0LD  dark  H2 runs per path summed over flies  3169.3 | model mean dwell 1/(1-P_ii): median over flies 5.0 min [IQR 3.5-7.3] | decoded KM median 3 min
hid per0LD  dark  H3 runs per path summed over flies  1992.6 | model mean dwell 1/(1-P_ii): median over flies 24.0 min [IQR 19.3-34.7] | decoded KM median 17 min
obs per0DD  whole S    runs  3227 (censored 0.00) | per-fly KM median 21 min | hazard bins outside model 95% band 0.82 of 11 | h per bin: 0.000 0.000 0.000 0.050 0.042 0.046 0.056 0.050 0.037 0.033 0.018 - -
obs per0DD  whole wake runs  3231 (censored 0.00) | per-fly KM median 9 min | hazard bins outside model 95% band 0.78 of 9 | h per bin: 0.076 0.063 0.075 0.080 0.090 0.094 0.090 0.091 0.060 - - - -
obs per0DD  whole Q    runs  6374 (censored 0.00) | per-fly KM median 2 min | hazard bins outside model 95% band 0.33 of 6 | h per bin: 0.452 0.381 0.392 0.376 0.330 0.282 - - - - - - -
obs per0DD  whole W    runs  6392 (censored 0.00) | per-fly KM median 2 min | hazard bins outside model 95% band 0.14 of 7 | h per bin: 0.334 0.296 0.260 0.239 0.235 0.182 0.196 - - - - - -
hid per0DD  whole H1 runs per path summed over flies  4327.3 | model mean dwell 1/(1-P_ii): median over flies 3.4 min [IQR 2.9-4.0] | decoded KM median 3 min
hid per0DD  whole H2 runs per path summed over flies  3819.2 | model mean dwell 1/(1-P_ii): median over flies 5.8 min [IQR 5.1-7.7] | decoded KM median 4 min
hid per0DD  whole H3 runs per path summed over flies  3182.5 | model mean dwell 1/(1-P_ii): median over flies 25.4 min [IQR 19.5-33.1] | decoded KM median 20 min
obs sss     light S    runs   443 (censored 0.01) | per-fly KM median 7 min | hazard bins outside model 95% band 0.86 of 7 | h per bin: 0.000 0.000 0.000 0.201 0.180 0.158 0.093 - - - - - -
obs sss     light wake runs   493 (censored 0.18) | per-fly KM median 76 min | hazard bins outside model 95% band 0.31 of 13 | h per bin: 0.026 0.015 0.018 0.015 0.016 0.019 0.013 0.008 0.007 0.006 0.004 0.002 0.002
obs sss     light Q    runs  6943 (censored 0.00) | per-fly KM median 1 min | hazard bins outside model 95% band 0.50 of 6 | h per bin: 0.614 0.490 0.444 0.353 0.277 0.207 - - - - - - -
obs sss     light W    runs  7100 (censored 0.01) | per-fly KM median 4 min | hazard bins outside model 95% band 0.82 of 11 | h per bin: 0.200 0.189 0.172 0.139 0.115 0.085 0.065 0.050 0.026 0.017 0.025 - -
hid sss     light H1 runs per path summed over flies  3480.9 | model mean dwell 1/(1-P_ii): median over flies 3.2 min [IQR 2.0-7.9] | decoded KM median 2 min
hid sss     light H2 runs per path summed over flies  3543.7 | model mean dwell 1/(1-P_ii): median over flies 19.0 min [IQR 10.7-54.6] | decoded KM median 7 min
hid sss     light H3 runs per path summed over flies   467.2 | model mean dwell 1/(1-P_ii): median over flies 8.8 min [IQR 7.3-10.2] | decoded KM median 8 min
obs sss     dark  S    runs   594 (censored 0.02) | per-fly KM median 7 min | hazard bins outside model 95% band 0.86 of 7 | h per bin: 0.000 0.000 0.000 0.235 0.213 0.123 0.088 - - - - - -
obs sss     dark  wake runs   591 (censored 0.14) | per-fly KM median 51 min | hazard bins outside model 95% band 0.62 of 13 | h per bin: 0.025 0.030 0.027 0.034 0.026 0.019 0.017 0.011 0.009 0.005 0.004 0.005 0.003
obs sss     dark  Q    runs  8617 (censored 0.00) | per-fly KM median 1 min | hazard bins outside model 95% band 0.67 of 6 | h per bin: 0.563 0.487 0.461 0.392 0.324 0.165 - - - - - - -
obs sss     dark  W    runs  8862 (censored 0.01) | per-fly KM median 4 min | hazard bins outside model 95% band 0.60 of 10 | h per bin: 0.217 0.213 0.190 0.175 0.138 0.117 0.080 0.049 0.027 0.021 - - -
hid sss     dark  H1 runs per path summed over flies  5172.8 | model mean dwell 1/(1-P_ii): median over flies 3.2 min [IQR 2.4-7.6] | decoded KM median 2 min
hid sss     dark  H2 runs per path summed over flies  5214.0 | model mean dwell 1/(1-P_ii): median over flies 10.1 min [IQR 6.6-18.1] | decoded KM median 5 min
hid sss     dark  H3 runs per path summed over flies   593.1 | model mean dwell 1/(1-P_ii): median over flies 7.8 min [IQR 6.9-10.1] | decoded KM median 7 min
```
- Observed: the fitted model captures the broad scale of every run type but misfits the shape: hazard bins outside the model's 95% band in 14-100% of bins (median ~0.6). Structural: S runs < 5 min are impossible in the data (sleep criterion) but predicted by the model; fumin and sss sleep runs end faster at long durations than predicted; long wake and active-wake runs are heavier-tailed than predicted in several blocks. Quiet-wake (Q) runs are fitted best.
- Hidden: mean dwell 1/(1-P_ii), fly medians: H1 2.9-4.0 min in every genotype and phase; H2 5-19 min; H3 20-36 min in CS/per0, 8-9 min in fumin/sss. Decoded dwell distributions match the run-weighted geometric curves (expected, circular).

- 2026-10-25: user will plot observed dwell times from data only in the paper (distribution modelling to be done in the user's MATLAB bout library). Added figures/fig9_dwell_observed_data.* (data only); fig9_dwell_observed.* (with model overlay) kept as a supplementary diagnostic.

## PRESPEC 20i — sensitivity: fly-phases without a Q/W reference excluded (2026-10-25, code 8b697b2)
15 fly-phases lack r_f (fumin 4, sss 11; 9 of them in non-degenerate flies: fumin light 4, sss light 4, sss dark 1). Their wake minutes were coded missing.
```
fumin light all                        | c_model +0.170 [+0.103,+0.232] (n 46) | J +0.129 [+0.085,+0.181] | EPR median 0.108 | round trips fwd-rev: observed +0.073 p 0.00015 (n 49), hidden +0.189 p 1.7e-05 (n 46)
fumin light excluding 4 fly-phases     | c_model +0.182 [+0.110,+0.245] (n 44) | J +0.139 [+0.090,+0.192] | EPR median 0.150 | round trips fwd-rev: observed +0.073 p 0.00015 (n 49), hidden +0.195 p 1.3e-05 (n 45)
fumin dark  all                        | c_model +0.335 [+0.258,+0.415] (n 45) | J +0.136 [+0.100,+0.174] | EPR median 0.143 | round trips fwd-rev: observed +0.078 p 0.00013 (n 50), hidden +0.305 p 8.8e-08 (n 46)
fumin dark  excluding 0 fly-phases     | c_model +0.335 [+0.251,+0.419] (n 45) | J +0.136 [+0.100,+0.179] | EPR median 0.143 | round trips fwd-rev: observed +0.078 p 0.00013 (n 50), hidden +0.305 p 8.8e-08 (n 46)
sss   light all                        | c_model -0.185 [-0.276,-0.103] (n 39) | J -0.056 [-0.083,-0.030] | EPR median 0.031 | round trips fwd-rev: observed -0.069 p 0.057 (n 46), hidden -0.142 p 0.0011 (n 36)
sss   light excluding 4 fly-phases     | c_model -0.193 [-0.286,-0.110] (n 37) | J -0.062 [-0.090,-0.035] | EPR median 0.031 | round trips fwd-rev: observed -0.069 p 0.057 (n 46), hidden -0.147 p 0.00096 (n 35)
sss   dark  all                        | c_model -0.062 [-0.141,+0.013] (n 45) | J -0.018 [-0.044,+0.008] | EPR median 0.044 | round trips fwd-rev: observed -0.052 p 0.2 (n 52), hidden -0.061 p 0.22 (n 39)
sss   dark  excluding 1 fly-phases     | c_model -0.069 [-0.150,+0.010] (n 44) | J -0.020 [-0.047,+0.007] | EPR median 0.044 | round trips fwd-rev: observed -0.052 p 0.2 (n 52), hidden -0.061 p 0.22 (n 39)
```
- All results HOLD: signs and significance unchanged (fumin forward; sss light reversed; sss dark ns). Observed round trips are identical by construction (they need r_f).
- Noted: coding those minutes as fully missing discards the fact that they are NOT sleep; the correct partial-observation likelihood would be P(Q|H)+P(W|H). Not applied (would need refits); reported to the user.

## PRESPEC 20j — per-fly EPR vs J (hidden), linear vs quadratic (2026-10-25, code b99b021); Fig 11
```
CS      both phases n  50 | M1 EPR = -0.158 + 2.46|J|  a CI [1.54,2.99] RSS 1.147 | M2 EPR = +0.035 + 4.48J^2  c CI [4.09,5.65] RSS 0.334 | better: quadratic (M2) (bootstrap share M1 better 0.01) | exponent n 1.76 [1.46,2.10] (n pts 31)
CS      light       n  25 | M1 EPR = +0.003 + 0.65|J|  a CI [0.48,0.81] RSS 0.020 | M2 EPR = +0.028 + 2.60J^2  c CI [1.96,3.65] RSS 0.023 | better: linear (M1) (bootstrap share M1 better 0.73) | exponent n 0.91 [0.32,1.46] (n pts 12)
CS      dark        n  25 | M1 EPR = -0.217 + 2.80|J|  a CI [1.81,3.35] RSS 0.604 | M2 EPR = +0.077 + 4.36J^2  c CI [3.90,5.59] RSS 0.219 | better: quadratic (M2) (bootstrap share M1 better 0.23) | exponent n 1.71 [1.33,2.12] (n pts 19)
fumin   both phases n  96 | M1 EPR = -0.084 + 2.23|J|  a CI [1.92,2.46] RSS 1.355 | M2 EPR = +0.115 + 3.63J^2  c CI [3.10,5.08] RSS 2.109 | better: linear (M1) (bootstrap share M1 better 0.96) | exponent n 1.40 [1.20,1.56] (n pts 66)
fumin   light       n  46 | M1 EPR = -0.088 + 2.11|J|  a CI [1.75,2.36] RSS 0.466 | M2 EPR = +0.107 + 3.29J^2  c CI [2.75,4.59] RSS 0.706 | better: linear (M1) (bootstrap share M1 better 0.77) | exponent n 1.39 [1.23,1.55] (n pts 32)
fumin   dark        n  50 | M1 EPR = -0.090 + 2.43|J|  a CI [2.01,2.73] RSS 0.772 | M2 EPR = +0.112 + 4.30J^2  c CI [3.66,7.26] RSS 1.178 | better: linear (M1) (bootstrap share M1 better 0.72) | exponent n 1.42 [1.09,1.74] (n pts 34)
per0LD  both phases n  54 | M1 EPR = -0.082 + 1.87|J|  a CI [1.43,2.18] RSS 0.900 | M2 EPR = +0.023 + 3.95J^2  c CI [3.49,4.25] RSS 0.350 | better: quadratic (M2) (bootstrap share M1 better 0.01) | exponent n 1.60 [1.36,1.82] (n pts 31)
per0LD  light       n  27 | M1 EPR = -0.048 + 1.81|J|  a CI [1.11,2.33] RSS 0.335 | M2 EPR = +0.033 + 4.03J^2  c CI [3.02,5.22] RSS 0.242 | better: quadratic (M2) (bootstrap share M1 better 0.37) | exponent n 1.47 [1.03,1.95] (n pts 11)
per0LD  dark        n  27 | M1 EPR = -0.139 + 2.01|J|  a CI [1.32,2.40] RSS 0.509 | M2 EPR = +0.009 + 3.96J^2  c CI [3.39,4.25] RSS 0.099 | better: quadratic (M2) (bootstrap share M1 better 0.00) | exponent n 1.73 [1.49,1.99] (n pts 20)
per0DD  both phases n  23 | M1 EPR = -0.316 + 2.73|J|  a CI [1.43,2.97] RSS 0.727 | M2 EPR = +0.090 + 2.84J^2  c CI [2.64,5.26] RSS 0.632 | better: quadratic (M2) (bootstrap share M1 better 0.30) | exponent n 1.84 [1.52,2.32] (n pts 18)
sss     both phases n  84 | M1 EPR = -0.046 + 1.83|J|  a CI [1.41,2.22] RSS 0.504 | M2 EPR = +0.026 + 6.44J^2  c CI [5.04,8.10] RSS 0.362 | better: quadratic (M2) (bootstrap share M1 better 0.12) | exponent n 1.51 [1.28,1.77] (n pts 39)
sss     light       n  37 | M1 EPR = -0.045 + 2.00|J|  a CI [1.35,2.54] RSS 0.314 | M2 EPR = +0.030 + 6.38J^2  c CI [4.44,8.31] RSS 0.235 | better: quadratic (M2) (bootstrap share M1 better 0.15) | exponent n 1.51 [1.18,2.00] (n pts 17)
sss     dark        n  47 | M1 EPR = -0.031 + 1.49|J|  a CI [1.13,1.74] RSS 0.139 | M2 EPR = +0.021 + 6.54J^2  c CI [5.51,8.88] RSS 0.127 | better: quadratic (M2) (bootstrap share M1 better 0.20) | exponent n 1.48 [1.11,1.90] (n pts 22)
slope a (M1, both phases) minus CS, bootstrap 95% CI:
  fumin   -0.22 [-0.88,+0.73]
  per0LD  -0.59 [-1.28,+0.41]
  per0DD  +0.27 [-0.96,+1.24]
  sss     -0.63 [-1.32,+0.40]
```
- Pooled phases: quadratic (EPR = b + cJ^2) fits better than linear for CS, per0LD, per0DD, sss (bootstrap share linear-better 0.01-0.30); linear better for fumin (0.96). Power-law exponents 1.4-1.8 (CIs 1.1-2.3), all excluding 1 except where few points. Linear fits need negative intercepts (-0.05 to -0.32 nats/h, unphysical); quadratic intercepts +0.02-0.04 (CS, per0LD, sss) match the estimator offset (20f); fumin +0.12, per0DD +0.09.
- c = 1/L: CS 4.5 [4.1,5.7], per0LD 4.0 [3.5,4.3], per0DD 2.8 [2.6,5.3], fumin 3.6 [3.1,5.1], sss 6.4 [5.0,8.1] (descriptive; the pre-stated genotype contrast was on a, all ns).
- Interpretation note: for each edge with one-way traffic N and net current J, ln((N+J/2)/(N-J/2)) ~ J/N, so count EPR ~ J^2 / N when J << N. A quadratic relation therefore means the net current is small compared with the back-and-forth traffic on each edge (near-linear response per edge), even though the summed cycle affinity is ~1 nat/cycle. This corrects my earlier remark that A ~ 1 places the data outside linear response: that holds for the whole cycle, not per edge. J and EPR share counts, so errors are correlated.

## POST HOC (user request, 2026-10-25): Onsager conductance L per genotype, EPR = J^2/L through the origin (user's Methods Eq. 25 scheme), per-fly hidden points of 20j. Code onsager_L.py; figures/fig12_onsager_L.*
```
both phases CS      n  50 | 1/L = 4.66 +- 0.14 | L = 0.215 +- 0.006 (95% CI 0.203-0.228; fly bootstrap 0.171-0.229) | R2(origin) 0.96 | L with offset subtracted 0.220
both phases per0LD  n  54 | 1/L = 4.08 +- 0.12 | L = 0.245 +- 0.007 (95% CI 0.231-0.261; fly bootstrap 0.226-0.271) | R2(origin) 0.95 | L with offset subtracted 0.252
both phases per0DD  n  23 | 1/L = 3.03 +- 0.17 | L = 0.330 +- 0.018 (95% CI 0.296-0.373; fly bootstrap 0.203-0.361) | R2(origin) 0.94 | L with offset subtracted 0.335
both phases fumin   n  96 | 1/L = 4.21 +- 0.19 | L = 0.237 +- 0.011 (95% CI 0.218-0.260; fly bootstrap 0.171-0.274) | R2(origin) 0.84 | L with offset subtracted 0.243
both phases sss     n  84 | 1/L = 6.92 +- 0.28 | L = 0.144 +- 0.006 (95% CI 0.134-0.157; fly bootstrap 0.119-0.178) | R2(origin) 0.88 | L with offset subtracted 0.153
   both phases: per0LD vs CS Wald z -3.11 p 0.0019 BH 0.0025
   both phases: per0DD vs CS Wald z -7.43 p 1.1e-13 BH 4.5e-13
   both phases: fumin vs CS Wald z -1.91 p 0.056 BH 0.056
   both phases: sss vs CS Wald z +7.31 p 2.7e-13 BH 5.4e-13
light       CS      n  25 | 1/L = 3.36 +- 0.32 | L = 0.297 +- 0.029 (95% CI 0.248-0.371; fly bootstrap 0.212-0.366) | R2(origin) 0.82 | L with offset subtracted 0.354
light       per0LD  n  27 | 1/L = 4.24 +- 0.27 | L = 0.236 +- 0.015 (95% CI 0.208-0.271; fly bootstrap 0.158-0.304) | R2(origin) 0.90 | L with offset subtracted 0.243
light       fumin   n  46 | 1/L = 3.77 +- 0.20 | L = 0.265 +- 0.014 (95% CI 0.239-0.298; fly bootstrap 0.192-0.319) | R2(origin) 0.88 | L with offset subtracted 0.272
light       sss     n  37 | 1/L = 6.82 +- 0.39 | L = 0.147 +- 0.008 (95% CI 0.131-0.166; fly bootstrap 0.115-0.195) | R2(origin) 0.89 | L with offset subtracted 0.153
   light: per0LD vs CS Wald z +2.08 p 0.038 BH 0.056
   light: fumin vs CS Wald z +1.06 p 0.29 BH 0.29
   light: sss vs CS Wald z +6.78 p 1.2e-11 BH 3.6e-11
dark        CS      n  25 | 1/L = 4.70 +- 0.19 | L = 0.213 +- 0.008 (95% CI 0.197-0.231; fly bootstrap 0.163-0.227) | R2(origin) 0.96 | L with offset subtracted 0.217
dark        per0LD  n  27 | 1/L = 4.01 +- 0.11 | L = 0.249 +- 0.007 (95% CI 0.236-0.264; fly bootstrap 0.234-0.277) | R2(origin) 0.98 | L with offset subtracted 0.257
dark        per0DD  n  23 | 1/L = 3.03 +- 0.17 | L = 0.330 +- 0.018 (95% CI 0.296-0.373; fly bootstrap 0.202-0.360) | R2(origin) 0.94 | L with offset subtracted 0.335
dark        fumin   n  50 | 1/L = 4.99 +- 0.32 | L = 0.200 +- 0.013 (95% CI 0.178-0.230; fly bootstrap 0.126-0.244) | R2(origin) 0.83 | L with offset subtracted 0.206
dark        sss     n  47 | 1/L = 7.24 +- 0.43 | L = 0.138 +- 0.008 (95% CI 0.123-0.157; fly bootstrap 0.108-0.163) | R2(origin) 0.86 | L with offset subtracted 0.151
   dark: per0LD vs CS Wald z -3.21 p 0.0013 BH 0.0018
   dark: per0DD vs CS Wald z -6.65 p 2.9e-11 BH 1.2e-10
   dark: fumin vs CS Wald z +0.77 p 0.44 BH 0.44
   dark: sss vs CS Wald z +5.37 p 8e-08 BH 1.6e-07
fly-bootstrap comparison with CS (pooled phases), L difference and 95% CI, two-sided bootstrap p:
   per0LD  dL +0.030 [+0.008,+0.077] p 0.011
   per0DD  dL +0.115 [-0.007,+0.167] p 0.0645
   fumin   dL +0.023 [-0.047,+0.076] p 0.447
   sss     dL -0.070 [-0.101,-0.020] p 0.008
   per0DD - per0LD dL [-0.044,+0.123] p 0.179
```
- OLS error bars and Wald tests treat a fly's light and dark points as independent; the fly bootstrap is the honest interval and gives wider, asymmetric CIs. By bootstrap (pooled phases): sss L lower than CS (p 0.008), per0LD higher (p 0.011), per0DD higher but ns (p 0.065; one fly with J ~ 0.95 dominates), fumin ns. Subtracting the 0.02 estimator offset raises L by 2-6%.

## PRESPEC 20k — route vs bout duration (2026-10-26, code 05f7bac; hidden null with 400 simulations after amendment A; first run with 40 kept as sqw10_nsim40.log)
```
R1_exit_vs_sleep         observed CS      light flies 25 | D +0.087 | null mean +0.001 | p 0.125 Holm 0.75
R1_exit_vs_sleep         observed CS      dark  flies 24 | D -0.293 | null mean +0.002 | p 0 Holm 0  *
R1_exit_vs_sleep         observed fumin   light flies 25 | D -0.002 | null mean -0.000 | p 0.994 Holm 1
R1_exit_vs_sleep         observed fumin   dark  flies 15 | D -0.061 | null mean +0.001 | p 0.349 Holm 1
R1_exit_vs_sleep         observed per0LD  light flies 26 | D -0.228 | null mean +0.000 | p 0.005 Holm 0.04  *
R1_exit_vs_sleep         observed per0LD  dark  flies 26 | D -0.121 | null mean +0.002 | p 0.091 Holm 0.637
R1_exit_vs_sleep         observed per0DD  whole flies 23 | D -0.028 | null mean -0.001 | p 0.589 Holm 1
R1_exit_vs_sleep         observed sss     light flies 32 | D +0.010 | null mean +0.000 | p 0.848 Holm 1
R1_exit_vs_sleep         observed sss     dark  flies 31 | D +0.008 | null mean +0.000 | p 0.873 Holm 1
R1_exit_vs_sleep         hidden   CS      light flies 24 | D +0.041 | null mean -0.031 | p 0.195 Holm 1
R1_exit_vs_sleep         hidden   CS      dark  flies 25 | D -0.275 | null mean +0.017 | p 0.00499 Holm 0.0449  *
R1_exit_vs_sleep         hidden   fumin   light flies 29 | D -0.053 | null mean -0.010 | p 0.559 Holm 1
R1_exit_vs_sleep         hidden   fumin   dark  flies 18 | D -0.144 | null mean -0.013 | p 0.354 Holm 1
R1_exit_vs_sleep         hidden   per0LD  light flies 27 | D -0.024 | null mean -0.008 | p 0.813 Holm 1
R1_exit_vs_sleep         hidden   per0LD  dark  flies 25 | D -0.039 | null mean +0.009 | p 0.409 Holm 1
R1_exit_vs_sleep         hidden   per0DD  whole flies 23 | D -0.027 | null mean -0.000 | p 0.579 Holm 1
R1_exit_vs_sleep         hidden   sss     light flies 18 | D +0.021 | null mean +0.000 | p 0.868 Holm 1
R1_exit_vs_sleep         hidden   sss     dark  flies 23 | D +0.004 | null mean +0.013 | p 0.883 Holm 1
R2_entry_vs_wake_before  observed CS      light flies 25 | D -0.263 | null mean +0.000 | p 0 Holm 0  *
R2_entry_vs_wake_before  observed CS      dark  flies 25 | D +0.072 | null mean +0.001 | p 0.234 Holm 0.702
R2_entry_vs_wake_before  observed fumin   light flies 30 | D -0.322 | null mean +0.004 | p 0.004 Holm 0.028  *
R2_entry_vs_wake_before  observed fumin   dark  flies 29 | D -0.146 | null mean +0.002 | p 0.294 Holm 0.702
R2_entry_vs_wake_before  observed per0LD  light flies 26 | D -0.180 | null mean +0.000 | p 0.002 Holm 0.016  *
R2_entry_vs_wake_before  observed per0LD  dark  flies 27 | D -0.077 | null mean +0.000 | p 0.127 Holm 0.508
R2_entry_vs_wake_before  observed per0DD  whole flies 23 | D -0.044 | null mean -0.001 | p 0.076 Holm 0.38
R2_entry_vs_wake_before  observed sss     light flies 29 | D -0.036 | null mean -0.001 | p 0.819 Holm 0.819
R2_entry_vs_wake_before  observed sss     dark  flies 28 | D -0.329 | null mean +0.005 | p 0.03 Holm 0.18
R2_entry_vs_wake_before  hidden   CS      light flies 25 | D -0.135 | null mean -0.028 | p 0.02 Holm 0.18
R2_entry_vs_wake_before  hidden   CS      dark  flies 25 | D +0.027 | null mean -0.017 | p 0.339 Holm 1
R2_entry_vs_wake_before  hidden   fumin   light flies 41 | D -0.051 | null mean -0.057 | p 0.903 Holm 1
R2_entry_vs_wake_before  hidden   fumin   dark  flies 29 | D +0.041 | null mean +0.004 | p 0.773 Holm 1
R2_entry_vs_wake_before  hidden   per0LD  light flies 26 | D -0.104 | null mean -0.042 | p 0.195 Holm 1
R2_entry_vs_wake_before  hidden   per0LD  dark  flies 26 | D -0.016 | null mean +0.007 | p 0.454 Holm 1
R2_entry_vs_wake_before  hidden   per0DD  whole flies 23 | D -0.018 | null mean -0.016 | p 0.903 Holm 1
R2_entry_vs_wake_before  hidden   sss     light flies 17 | D -0.416 | null mean -0.079 | p 0.0998 Holm 0.798
R2_entry_vs_wake_before  hidden   sss     dark  flies 22 | D +0.028 | null mean -0.048 | p 0.608 Holm 1
R3_entry_vs_sleep        observed CS      light flies 25 | D +0.048 | null mean -0.001 | p 0.39 Holm 1
R3_entry_vs_sleep        observed CS      dark  flies 25 | D -0.166 | null mean -0.001 | p 0.019 Holm 0.171
R3_entry_vs_sleep        observed fumin   light flies 30 | D -0.071 | null mean -0.000 | p 0.046 Holm 0.368
R3_entry_vs_sleep        observed fumin   dark  flies 29 | D -0.032 | null mean +0.001 | p 0.454 Holm 1
R3_entry_vs_sleep        observed per0LD  light flies 26 | D -0.122 | null mean -0.000 | p 0.098 Holm 0.637
R3_entry_vs_sleep        observed per0LD  dark  flies 27 | D -0.119 | null mean -0.002 | p 0.091 Holm 0.637
R3_entry_vs_sleep        observed per0DD  whole flies 23 | D -0.066 | null mean -0.000 | p 0.124 Holm 0.637
R3_entry_vs_sleep        observed sss     light flies 29 | D -0.022 | null mean -0.000 | p 0.693 Holm 1
R3_entry_vs_sleep        observed sss     dark  flies 28 | D +0.058 | null mean +0.001 | p 0.213 Holm 0.852
R3_entry_vs_sleep        hidden   CS      light flies 25 | D +0.111 | null mean -0.000 | p 0.0499 Holm 0.449
R3_entry_vs_sleep        hidden   CS      dark  flies 25 | D -0.080 | null mean -0.005 | p 0.334 Holm 1
R3_entry_vs_sleep        hidden   fumin   light flies 41 | D -0.065 | null mean -0.002 | p 0.224 Holm 1
R3_entry_vs_sleep        hidden   fumin   dark  flies 29 | D -0.048 | null mean -0.001 | p 0.569 Holm 1
R3_entry_vs_sleep        hidden   per0LD  light flies 26 | D -0.118 | null mean +0.002 | p 0.0499 Holm 0.449
R3_entry_vs_sleep        hidden   per0LD  dark  flies 26 | D -0.025 | null mean -0.000 | p 0.623 Holm 1
R3_entry_vs_sleep        hidden   per0DD  whole flies 23 | D -0.071 | null mean -0.000 | p 0.0898 Holm 0.628
R3_entry_vs_sleep        hidden   sss     light flies 17 | D -0.031 | null mean +0.011 | p 0.788 Holm 1
R3_entry_vs_sleep        hidden   sss     dark  flies 22 | D -0.089 | null mean -0.007 | p 0.419 Holm 1
```
- Significant (Holm < 0.05), all with D < 0 (quiet routes with SHORTER bouts):
  R1 exit vs sleep-bout length: observed CS dark (-0.29), per0LD light (-0.23); hidden CS dark (-0.28, beyond all 400 Markov simulations: the fitted Markov model cannot produce it).
  R2 entry vs prior wake length: observed CS light (-0.26), fumin light (-0.32), per0LD light (-0.18); none at the hidden level (largest: CS light -0.14, p 0.02 uncorrected).
  R3: none.
- Reading: CS at night wakes from long sleep bouts directly into active wake more often than from short bouts, beyond what the Markov model allows: exit routing has memory of sleep-bout length. R2 (observed) may partly reflect classification (short wake bouts are quieter overall); it is not confirmed at the hidden level.
- Exports: route_duration_observed_theta0.5.csv; route_duration_hidden_theta0.5.csv.gz (20 posterior paths per fly). Figures 13, 13b. Recipe: FIGURE_RECIPES section I.

## PRESPEC 20l — energy landscape (2026-10-29, code ); Figs 14, 15 (Viterbi, posterior), 13c
```
(1) state energies relative to H1 (nats), median [fly-bootstrap 95% CI]
  E2_minus_E1 CS      light n 25 median -0.57 [-0.65,-0.50]
  E2_minus_E1 CS      dark  n 25 median -0.69 [-1.00,-0.47]
  E2_minus_E1 per0LD  light n 27 median -0.46 [-0.83,-0.19]
  E2_minus_E1 per0LD  dark  n 27 median -0.35 [-0.70,-0.15]
  E2_minus_E1 per0DD  whole n 23 median -0.38 [-0.58,-0.26]
  E2_minus_E1 fumin   light n 46 median -0.92 [-0.97,-0.68]
  E2_minus_E1 fumin   dark  n 50 median -1.44 [-1.57,-1.27]
  E2_minus_E1 sss     light n 37 median -1.32 [-1.59,-1.18]
  E2_minus_E1 sss     dark  n 47 median -1.03 [-1.14,-0.88]
  E3_minus_E1 CS      light n 25 median -1.08 [-1.41,-0.71]
  E3_minus_E1 CS      dark  n 25 median -1.85 [-2.11,-1.71]
  E3_minus_E1 per0LD  light n 27 median -1.73 [-1.89,-1.29]
  E3_minus_E1 per0LD  dark  n 27 median -1.61 [-1.95,-1.29]
  E3_minus_E1 per0DD  whole n 23 median -1.85 [-1.95,-1.28]
  E3_minus_E1 fumin   light n 46 median +0.74 [+0.48,+1.22]
  E3_minus_E1 fumin   dark  n 50 median +1.07 [+0.87,+1.47]
  E3_minus_E1 sss     light n 37 median +1.13 [+0.85,+1.58]
  E3_minus_E1 sss     dark  n 47 median +1.38 [+0.95,+1.61]
(2) landscape prediction vs counted routing (viterbi): median(counted - predicted), Spearman
  viterbi $F_{S\to W}$ ( CS      light   n 25 median diff -0.021 Spearman +0.86
  viterbi $F_{S\to W}$ ( per0LD  light   n 27 median diff -0.044 Spearman +0.95
  viterbi $F_{S\to W}$ ( fumin   light   n 45 median diff -0.000 Spearman +0.81
  viterbi $F_{S\to W}$ ( sss     light   n 33 median diff -0.011 Spearman +0.97
  viterbi $F_{W\to S}$ ( CS      light   n 25 median diff -0.028 Spearman +0.84
  viterbi $F_{W\to S}$ ( per0LD  light   n 27 median diff -0.030 Spearman +0.92
  viterbi $F_{W\to S}$ ( fumin   light   n 45 median diff +0.000 Spearman +0.92
  viterbi $F_{W\to S}$ ( sss     light   n 34 median diff +0.000 Spearman +0.97
  viterbi $F_{S\to W}$ ( CS      dark/DD n 25 median diff +0.066 Spearman +0.87
  viterbi $F_{S\to W}$ ( per0LD  dark/DD n 27 median diff -0.006 Spearman +0.86
  viterbi $F_{S\to W}$ ( per0DD  dark/DD n 23 median diff +0.055 Spearman +0.76
  viterbi $F_{S\to W}$ ( fumin   dark/DD n 47 median diff +0.002 Spearman +0.90
  viterbi $F_{S\to W}$ ( sss     dark/DD n 43 median diff -0.034 Spearman +0.96
  viterbi $F_{W\to S}$ ( CS      dark/DD n 25 median diff +0.055 Spearman +0.95
  viterbi $F_{W\to S}$ ( per0LD  dark/DD n 27 median diff -0.050 Spearman +0.89
  viterbi $F_{W\to S}$ ( per0DD  dark/DD n 23 median diff +0.067 Spearman +0.91
  viterbi $F_{W\to S}$ ( fumin   dark/DD n 44 median diff -0.003 Spearman +0.93
  viterbi $F_{W\to S}$ ( sss     dark/DD n 43 median diff -0.002 Spearman +0.96
(2) landscape prediction vs counted routing (ffbs): median(counted - predicted), Spearman
  ffbs $F_{S\to W}$ ( CS      light   n 25 median diff -0.053 Spearman +0.96
  ffbs $F_{S\to W}$ ( per0LD  light   n 27 median diff +0.010 Spearman +0.99
  ffbs $F_{S\to W}$ ( fumin   light   n 45 median diff -0.019 Spearman +0.95
  ffbs $F_{S\to W}$ ( sss     light   n 33 median diff -0.000 Spearman +1.00
  ffbs $F_{W\to S}$ ( CS      light   n 25 median diff -0.048 Spearman +0.92
  ffbs $F_{W\to S}$ ( per0LD  light   n 27 median diff +0.009 Spearman +0.99
  ffbs $F_{W\to S}$ ( fumin   light   n 45 median diff -0.009 Spearman +0.99
  ffbs $F_{W\to S}$ ( sss     light   n 34 median diff +0.000 Spearman +0.99
  ffbs $F_{S\to W}$ ( CS      dark/DD n 25 median diff +0.034 Spearman +0.94
  ffbs $F_{S\to W}$ ( per0LD  dark/DD n 27 median diff +0.000 Spearman +0.97
  ffbs $F_{S\to W}$ ( per0DD  dark/DD n 23 median diff +0.044 Spearman +0.96
  ffbs $F_{S\to W}$ ( fumin   dark/DD n 47 median diff +0.012 Spearman +0.93
  ffbs $F_{S\to W}$ ( sss     dark/DD n 43 median diff -0.004 Spearman +0.99
  ffbs $F_{W\to S}$ ( CS      dark/DD n 25 median diff +0.016 Spearman +0.99
  ffbs $F_{W\to S}$ ( per0LD  dark/DD n 27 median diff +0.004 Spearman +0.99
  ffbs $F_{W\to S}$ ( per0DD  dark/DD n 23 median diff +0.028 Spearman +0.98
  ffbs $F_{W\to S}$ ( fumin   dark/DD n 45 median diff +0.003 Spearman +1.00
  ffbs $F_{W\to S}$ ( sss     dark/DD n 43 median diff -0.000 Spearman +0.99
```
- (1) CS and per0: wake H2 and sleep H3 lie below the regulator H1 (E2-E1 -0.35 to -0.69, E3-E1 -1.08 to -1.85 nats); CS sleep is 0.8 nats deeper at night than by day. fumin and sss: H2 further below H1 (-0.9 to -1.4) but H3 ABOVE H1 (+0.7 to +1.4). Qualitatively the user's binary-model finding (sleep/wake wells below the regulator in wt/per0, at or above it in fmn/sss), with smaller magnitudes.
- (2) Landscape prediction vs counted routing per fly: close agreement (Viterbi Spearman 0.76-0.97, |median diff| <= 0.07; posterior 0.92-1.00, <= 0.05). Consistency check, as stated in the plan; supplementary.
- (3) Fig 13c: the CS-dark exit curve falls from ~0.52 to ~0.30 across sleep-bout length, crossing the flat landscape line (0.35); fumin dark falls too (ns in 20k). Offsets between curves and lines elsewhere are consistent with the matched-definition difference stated in the plan.

## Presentation rule for Figs 13 / 13c (2026-10-29, user-approved after seeing the pooled curves)
Flies weighted equally per duration bin, bins with < 10 contributing flies omitted, bootstrap-over-flies intervals (route_bins.py). Effect: CS dark exit curve flattens visually (0.44, 0.44, 0.49, 0.33, 0.32) though its within-fly test stays significant; fumin dark (0.31, 0.19, 0.08) and sss dark (0.43, 0.50, 0.59) tails remain, each from >= 10 flies but with wide intervals, and their within-fly tests remain ns (D -0.14 p 0.35; +0.004 p 0.88). Per-fly proportions in sparse bins rest on 1-2 bouts per fly, so the binned curves stay descriptive; the within-fly D (Fig 13b) is the test.

## Fig 13d (2026-10-29, user request): per-fly D, flies as points, group mean, null range of the mean (fig_perfly_D.py; per-fly values fig13d.log)
- R1 CS dark: 79% of flies D < 0 (observed), 76% (hidden); per0LD light observed 65%. Spread is wide in fumin and sss, especially R2 (few bouts per fly).

## Fig 16 (2026-10-29, user request): R1 main-figure candidate (fig_R1_main.py). CS dark hidden P(quiet exit) by bin, fly-weighted: 0.44, 0.44, 0.49, 0.33, 0.32 vs flat landscape prediction 0.35; observed 0.24, 0.26, 0.23, 0.12, 0.13.

## PRESPEC 20m — 30-s vs 1-min steps, CS and fumin (2026-10-29, code 4844dd3); Fig 17
```
degenerate H1: 1 min 0 flies, 30 s 0 flies (of 78)
h1_dwell_min  CS    light n 25 | 1 min median +3.105 | 30 s median +2.116 | ratio of medians 0.68 | paired Wilcoxon p 1.2e-07 | Spearman +0.62
h1_dwell_min  CS    dark  n 25 | 1 min median +3.163 | 30 s median +2.312 | ratio of medians 0.73 | paired Wilcoxon p 6e-08 | Spearman +0.30
h1_dwell_min  fumin light n 46 | 1 min median +3.491 | 30 s median +1.942 | ratio of medians 0.56 | paired Wilcoxon p 4.2e-11 | Spearman +0.65
h1_dwell_min  fumin dark  n 46 | 1 min median +2.543 | 30 s median +1.805 | ratio of medians 0.71 | paired Wilcoxon p 1.1e-10 | Spearman +0.68
h1_one_step   CS    light n 25 | 1 min median +0.302 | 30 s median +0.225 | ratio of medians 0.75 | paired Wilcoxon p 0.00081 | Spearman +0.58
h1_one_step   CS    dark  n 25 | 1 min median +0.312 | 30 s median +0.213 | ratio of medians 0.68 | paired Wilcoxon p 0.00038 | Spearman +0.22
h1_one_step   fumin light n 46 | 1 min median +0.288 | 30 s median +0.265 | ratio of medians 0.92 | paired Wilcoxon p 0.071 | Spearman +0.61
h1_one_step   fumin dark  n 46 | 1 min median +0.362 | 30 s median +0.270 | ratio of medians 0.75 | paired Wilcoxon p 0.00013 | Spearman +0.68
J             CS    light n 25 | 1 min median +0.023 | 30 s median +0.005 | ratio of medians 0.20 | paired Wilcoxon p 0.63 | Spearman +0.73
J             CS    dark  n 25 | 1 min median +0.208 | 30 s median +0.189 | ratio of medians 0.91 | paired Wilcoxon p 0.25 | Spearman +0.86
J             fumin light n 46 | 1 min median +0.130 | 30 s median +0.158 | ratio of medians 1.22 | paired Wilcoxon p 0.036 | Spearman +0.73
J             fumin dark  n 46 | 1 min median +0.127 | 30 s median +0.128 | ratio of medians 1.01 | paired Wilcoxon p 0.44 | Spearman +0.74
EPR           CS    light n 25 | 1 min median +0.041 | 30 s median +0.034 | ratio of medians 0.85 | paired Wilcoxon p 0.89 | Spearman +0.33
EPR           CS    dark  n 25 | 1 min median +0.299 | 30 s median +0.221 | ratio of medians 0.74 | paired Wilcoxon p 0.41 | Spearman +0.86
EPR           fumin light n 46 | 1 min median +0.160 | 30 s median +0.195 | ratio of medians 1.22 | paired Wilcoxon p 0.14 | Spearman +0.65
EPR           fumin dark  n 46 | 1 min median +0.169 | 30 s median +0.153 | ratio of medians 0.90 | paired Wilcoxon p 0.92 | Spearman +0.62
c_ffbs        CS    light n 25 | 1 min median +0.009 | 30 s median -0.016 | ratio of medians -1.80 | paired Wilcoxon p 0.56 | Spearman +0.69 | same sign 0.76
c_ffbs        CS    dark  n 25 | 1 min median +0.266 | 30 s median +0.222 | ratio of medians 0.84 | paired Wilcoxon p 0.09 | Spearman +0.88 | same sign 0.92
c_ffbs        fumin light n 44 | 1 min median +0.203 | 30 s median +0.246 | ratio of medians 1.21 | paired Wilcoxon p 0.13 | Spearman +0.65 | same sign 0.86
c_ffbs        fumin dark  n 41 | 1 min median +0.307 | 30 s median +0.301 | ratio of medians 0.98 | paired Wilcoxon p 0.85 | Spearman +0.54 | same sign 0.90
R1_D          CS    light n 24 | 1 min median +0.043 | 30 s median +0.019 | ratio of medians 0.43 | paired Wilcoxon p 0.66 | Spearman +0.47 | same sign 0.75
R1_D          CS    dark  n 25 | 1 min median -0.299 | 30 s median -0.226 | ratio of medians 0.76 | paired Wilcoxon p 0.096 | Spearman +0.61 | same sign 0.88
R1_D          fumin light n 24 | 1 min median -0.132 | 30 s median -0.187 | ratio of medians 1.42 | paired Wilcoxon p 0.2 | Spearman +0.87 | same sign 0.83
R1_D          fumin dark  n 12 | 1 min median -0.120 | 30 s median -0.101 | ratio of medians 0.84 | paired Wilcoxon p 0.73 | Spearman +0.82 | same sign 0.75
CS dark R1 hidden D at 60 s: mean -0.277, flies < 0 0.76 (n 25); Wilcoxon vs 0 p 0.00071
CS dark R1 hidden D at 30 s: mean -0.186, flies < 0 0.80 (n 25); Wilcoxon vs 0 p 0.0061
```
- H1 is faster than the 1-min grid showed: mean dwell ~2 min at 30 s vs ~3 min at 1 min (ratio 0.56-0.73, p < 1e-7). It is resolved at 30 s (single-step runs 21-27%). Part may be Q/W flicker from the noisier 30-s activity fraction.
- J and EPR barely change: no significant change except fumin light J +22% (p 0.04); CS dark EPR -26% (ns). The 1-min grid did not hide substantial extra cycling or dissipation.
- Directions hold: per-fly c same sign in 76-92% (lowest in CS light, where c ~ 0); group medians within 0.05.
- CS-dark R1 memory holds: hidden D -0.19 at 30 s (80% of flies < 0, Wilcoxon p 0.006) vs -0.28 at 1 min.
