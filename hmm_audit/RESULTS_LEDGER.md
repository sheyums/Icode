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
