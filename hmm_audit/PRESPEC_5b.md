# PRESPEC 5b — light/dark split, added BEFORE any PRESPEC 5 output was viewed
Date 2026-10-19. pooled.py (PRESPEC 5) is running; its log has not been opened.
New information from the user: dt = 1 s, recording starts at ZT0, LD 12:12, windows ~4 h.
- PRESPEC 5's whole-record analysis remains PRIMARY and unchanged.
- Secondary: repeat P1 and P2 within each phase, LIGHT (ZT 0-12) and DARK (ZT 12-24),
  per genotype. No other split (no per-day, no per-window-of-day).
- Time of a bout = its START frame: frames before it in the fly's record, counting each empty
  column as 14400 frames (nominal 4 h, since its length is unknown). ZT = (t/3600) mod 24.
  Phase is assigned by bout START (pipeline convention).
- Bouts crossing a ZT0/ZT12 boundary carry out-of-phase time and are NOT censored (no
  censoring support): report, per genotype x phase x class, the fraction of bouts that cross.
- P2 pairs are assigned to the phase of their FIRST bout. Same permutation/bootstrap
  settings as PRESPEC 5 (2000 each, seed 5); P1 same K range, starts and CV folds.
- Light-vs-dark differences in A: bootstrap CI of the difference (same flies resampled for
  both phases, since phases share flies) excludes 0.
