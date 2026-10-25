# PRESPEC 20g — round-trip types around each sleep bout, observed and hidden (before running)
Date 2026-10-25. User request, intended for main figures. Descriptive in origin (the per-transition rates behind it
were seen), but the definitions and tests are fixed here before counting.

## Round trip = one sleep bout with its entry and exit
Four types, named by entry/exit: forward = quiet in, active out (Q->S->A; hidden H1->H3->H2); reverse = active in,
quiet out (A->S->Q; H2->H3->H1); direct = A->S->A (H2->H3->H2); quiet-both = Q->S->Q (H1->H3->H1).
Phase of a round trip = phase (light ZT 0-12 / dark ZT 12-24; per0DD whole) of the sleep-bout ONSET.
- Observed (1-s counts, PRESPEC 19 classification, bridging (b), theta 0.5): sleep bouts whose preceding AND following
  wake bouts each last >= 10 min. Entry is Q if the active fraction in seconds 61-300 before onset < theta r_f, exit
  is Q if the active fraction in seconds 61-300 after waking < theta r_f; r_f = fly reference of the onset's phase
  (sqw_ref.wake_reference). Flies need >= 3 round trips in the phase.
- Hidden (20b fits, theta 0.5): on each of 200 FFBS paths, every H3 run whose neighbours are both H1/H2 (not the
  record ends); entry = state before, exit = state after. Fractions per path, averaged over paths. Non-degenerate
  flies with >= 3 H3 runs per phase on average.
## Reported per genotype x phase
Mean fraction of each type over flies, 95% fly-bootstrap interval (2000, seed 20); fraction of flies whose most
common type is forward; forward - reverse per fly, Wilcoxon signed-rank vs 0 (Holm over the 9 blocks, separately for
observed and hidden). Most common type per block = largest mean fraction.
Sensitivity: theta 0.33 and 0.75 for both levels (hidden from the K = 3 theta fits).
## Caveat fixed now
The observed level depends on the Q threshold: theta moves every type's share; the forward-minus-reverse difference is
the quantity robust to theta (PRESPEC 20b theta check). Statements about which type is most common are reported with
their theta range.
