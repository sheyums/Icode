# PRESPEC 20l — energy landscape: state energies, landscape vs counted routing, landscape vs duration (before running)
Date 2026-10-29. User request; descriptive. 20b fits (theta 0.5), labelled H1, H2, H3; non-degenerate flies; fly-phases
without a Q/W reference or with no cycle excluded (as 20i/20j).
## (1) State energies relative to the regulator
Per fly x phase, from the phase matrix's stationary distribution pi (user's Methods Eqs. 5, 7):
E2 - E1 = ln(pi_1 / pi_2), E3 - E1 = ln(pi_1 / pi_3) (nats; E1 = 0). Plot: individual flies as points (jittered
horizontally), median as a short horizontal bar, no box; genotypes on x, light (open) and dark/DD (filled) side by
side; two panels (E2 - E1, E3 - E1). Report medians with fly-bootstrap 95% intervals.
## (2) Landscape prediction vs counted routing
x = landscape prediction per fly x phase: F_S->W = tanh[(E32 - E31c)/2] = (r - 1)/(r + 1), r = P31 b12 / P32,
b12 = P12/(P12 + P13); F_W->S likewise with P21 b13 / P23 (Methods Eqs. 10-16, from the phase matrix).
y = counted routing from the same fly-phase: (a) Viterbi path counts (user request), (b) FFBS posterior counts
(per-path Xn/Yn averaged, >= 3 routes per denominator), as a second figure.
Layout: 2 x 2 (columns light, dark/DD; rows F_S->W = Xn, F_W->S = Yn), all genotypes, identity line, axes -1..1.
Reported: per genotype x phase, median (counted - predicted) and Spearman correlation. Stated now: the fitted model is
estimated from the same data, so agreement on fly averages is expected; this is a consistency check.
## (3) Landscape prediction vs bout duration (Fig. 13 with landscape lines)
Hidden-level panels of Fig. 13 (R1 quiet exit vs sleep-bout length; R2/R3 quiet entry), with a horizontal line per
genotype x phase at the landscape prediction converted to a probability: P(quiet exit) = (F_S->W + 1)/2,
P(quiet entry) = (F_W->S + 1)/2, mean over flies. Stated now: the prediction is for all exits / entries, the counted
curves use the matched definition (H3 runs >= 5 min, gaps >= 10 min), so a constant offset between them can come
from the definitions; the informative feature is a SLOPE of the counted curve against a flat prediction.
