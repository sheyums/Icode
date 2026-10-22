# PRESPEC 17 — per-fly E2 effect sizes and the pre-sleep taper (before running)
Date 2026-10-22. Count data, all 5 groups. Definitions of sleep onset, waking and movement fraction as
PRESPEC 15. Blocks as PRESPEC 16 (light/dark halves; per0DD whole record only).
## A. Per-fly asymmetry
Wake bouts >= 20 min. Per fly x block: D = mean over bouts and 1-min bins (min 2-10) of (EXIT - APPROACH),
flies with >= 3 eligible bouts. Per fly: permutation p (2000 swaps of EXIT/APPROACH per bout).
Per group x block: median D, IQR, fraction of flies with D > 0, exact two-sided sign test across flies,
fraction of flies with permutation p < 0.05.
## B. Taper depth (per fly)
Wake bouts >= 20 min. Reference = movement fraction in the middle 5 min of the bout.
Depth = reference - movement fraction in minutes 2-3 before sleep onset; relative depth = depth / reference
(flies with reference > 0). Per group x block: median and IQR over flies; Wilcoxon signed-rank of depth vs 0.
## C. Relaxation times into and out of sleep (group level)
Wake bouts >= 40 min. Approach profile a(t), t = 2..20 min before onset; exit profile e(t), t = 2..20 min
after waking; each divided by that bout's mid-bout reference (middle 5 min), averaged over bouts within
fly, then over flies. Fit y(t) = 1 - A exp(-(t-1)/tau) by least squares (A in [-5,5], tau in [0.1, 100] min)
to each mean profile: tau_in (approach) and tau_out (exit). 1000 fly-bootstrap resamples for 95% intervals.
Time-reversal symmetry predicts tau_in = tau_out and A_in = A_out; reported: the difference with its
bootstrap interval. Groups x blocks with < 8 flies contributing are not fitted.
## Figures
(i) diagnostic report section (all numbers, all blocks); (ii) journal-style figures (vector PDF + PNG):
normalised approach/exit profiles per group; per-fly D distribution; tau_in vs tau_out with intervals.
