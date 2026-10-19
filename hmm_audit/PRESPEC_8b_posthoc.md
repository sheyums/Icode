# 8b — POST HOC robustness check of PRESPEC 8 (written AFTER seeing PRESPEC 8 results)
Concern: windows were cut at bout boundaries near each 4 h mark. If the cut rule is asymmetric in
time (e.g. first boundary AFTER the mark), the last bout of a window is length-biased and the first
is not. That is a time-asymmetric selection the per-window reversal null cannot absorb, so it can
produce D* > 0 from a reversible process.
Checks:
1. Window cut position: (cumulative frames at window end) - k*14400, distribution, per genotype.
2. Mean log-duration of the FIRST, INTERIOR and LAST bout of each window, per class and genotype.
3. Rerun PRESPEC 8 exactly, but using only pairs whose two bouts are both NOT the first or last bout
   of their window (edge-free pairs). Same bins, null and bootstrap.
Interpretation: PRESPEC 8 conclusions stand only for blocks that still pass in check 3.
