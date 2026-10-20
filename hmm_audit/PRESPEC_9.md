# PRESPEC 9 — two kinds of sleep bout within a fly, EDGE-FREE, at 4 h and 12 h blocks (before running)
Date 2026-10-19. Written AFTER PRESPEC 7 results (CS supported 10/94 under N2; per0, fumin not;
sss 4 blocks) and after 8b found the window-edge length bias. User confirmed the cut rule:
windows were extended so a bout already under way is not cut, so the LAST bout of each window is
length-selected. The 12 h blocks are a WIDER block chosen after seeing the 4 h result; they are
reported ALONGSIDE the 4 h result and never replace it (Icode rule 4).

## Bouts
- Sleep bouts as PRESPEC 5/7 (interior runs, xmin 300 s, dt 1 s, ZT from start frame with empty
  columns = 14400 frames), EXCLUDING the last run of every window column (the length-selected one).

## A. 4 h blocks (rerun of PRESPEC 7 on edge-free bouts)
- Block = fly x ZT bin [0,4)...[20,24), days pooled, n >= 20. LR test K=1 vs K=2 via em_fit.fit2g.
- Nulls, 100 parametric replicates each: N1 one rate; N2 one rate per day (PRIMARY null).
## B. 12 h blocks
- Block = fly x phase (light ZT0-12, dark ZT12-24), days pooled, n >= 20.
- Nulls, 100 replicates each: N1 one rate; N2 one rate per day; N3 one rate per (day x 4 h sub-bin)
  (PRIMARY null for B: absorbs within-phase time-of-day drift as well as day-to-day drift).
## Read-out (both A and B, per genotype)
- PRIMARY: Fisher's combined p over the genotype's blocks, using each block's primary-null p.
  Supported if Fisher p < 0.05. (Bootstrap p-values are discrete, min 1/101; noted, not corrected.)
- Also reported: count of blocks with p < 0.05 and its exact binomial p against 5%, for every null.
- For blocks with primary p < 0.05: median short-component weight and the two component means.
- A and B are both reported; if they disagree, A answers "within one fly at one time of day".
- Seeds: SeedSequence([9, A/B, genotype, fly, bin]).
