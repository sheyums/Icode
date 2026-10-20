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
