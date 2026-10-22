# PRESPEC 19 — signed cycle current in observable states (before running)
Date 2026-10-23. Requested by the user after Fig. 3 suggested sleepless sits above the diagonal. The +-5 min
window was chosen after inspection (PRESPEC 18), so all p-values here are DESCRIPTIVE; any claim of a reversed
sleepless current needs independent flies.
Data: count files, all 5 groups. Sleep = inactive run >= 300 s with bridging (b) (a 1 s active run bridged only
when both flanking inactive runs are >= 300 s; confirmed lab convention). Blocks as PRESPEC 16 (light/dark by
transition time; per0DD whole record).
## Observable states at each transition
- Each sleep ONSET is preceded by a wake bout; each WAKING is followed by one. Only wake bouts >= 10 min are
  used, so the two +-5 min windows of a bout do not overlap; their onset and waking are both classified.
- Transit activity = fraction of active seconds in seconds 61-300 before the onset (approach) or after the
  waking (exit). Fly reference r_f = median, over that fly's eligible wake bouts in the block, of activity in the
  bout's middle 5 min.
- Q (quiet wake) if transit activity < theta * r_f, else A (active wake). PRIMARY theta = 0.5; sensitivity
  theta = 0.33 and 0.75, all reported.
## Current
Forward cycle A -> Q -> S -> A. Per fly: n_QS, n_AS (onsets via Q / from A), n_SQ, n_SA (wakings into Q / A).
- Normalised current c = n_QS/(n_QS+n_AS) - n_SQ/(n_SQ+n_SA)  (= P(Q|onset) - P(Q|waking)), in [-1, 1].
- Rate J = [(n_QS - n_SQ) + (n_SA - n_AS)] / 2 per hour of the block's recording time (cycles per hour).
Flies with >= 3 eligible wake bouts.
## Statistics (descriptive)
Per group x block: mean c with 95% fly-bootstrap interval (2000), median, fraction of flies with c > 0, exact
two-sided sign test, Wilcoxon signed-rank vs 0; same for J. Contrasts: sleepless vs each other genotype in the
same block (Mann-Whitney on per-fly c). Seed 19.
## Figure
c (and J) per group x block, mean +- 95% interval, per-fly points, plus a cycle diagram; journal style.
