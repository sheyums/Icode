# PRESPEC 20w-B — clock-face cosine test (replaces 20w's summary statistic only; before any run)
Date 2026-11-01. Owner: cClaudeII (Plan B feasibility session), on cClaude's brief relayed by Sheyum. Branch
planB-feasibility. Results go to RESULTS_LEDGER_planB.md.

## Why this exists, stated openly
20w's pre-registered summary, the range (max - min) of the 12-bin J profile, has no power. Noise inflates the range
of noisy bins, and the bootstrap inflates it further: per0 DD's band [0.120, 0.355] lay entirely above its own
point estimate 0.092, so no genotype could pass. The 20w result STANDS as recorded (descriptive; no genotype
rhythmic by its rule). This statistic was chosen AFTER seeing that the range was insensitive. The reason is power,
not the direction of any effect. Data, bins and pooling are 20w's, unchanged. No re-binning (root CLAUDE.md rule 4).

## Data (fixed now)
- Primary family (Holm): CS LD, CS DD, per0 LD.
- Negative control: per0 DD.
- Descriptive only, outside the family: fumin, sss.
Each fly is its 20w record: sqw_ref.sqw_symbols(theta 0.5; LD light_dark=True, DD False), trimmed to whole 24-h days
(sqw_symbols' full_days). LD days are aligned at ZT0. CS DD is binned in clock time from the cut start. The
free-running smear (period spread across flies, up to ~2 h over the record) is a known dilution of any CS DD
amplitude.

## Statistic
Per genotype, 20w's pooling. Transitions between consecutive non-missing minutes are attributed to the 2-h clock bin
b = 0..11 of the arrival minute, pooled over flies and days. Bin b's pooled first-order chain gives J (cycles/h,
sqw_ref.model_energetics). First harmonic by least squares over the bins with a defined J:
J(b) = m + a cos(2 pi b / 12) + c sin(2 pi b / 12). Statistic A = sqrt(a^2 + c^2).
Phase atan2(c, a), reported as clock time of the peak (bin START convention, +1 h for the bin centre); descriptive.
Secondary, outside the Holm family: the same fit to the pooled c_obs profile (onset/waking events pooled per bin,
as 20w). EPR profile: figure only, no test (zero-count edges make it unstable bin by bin).

## Null (one per genotype)
B = 999 replicates; replicate k uses default_rng(20_900 + k) for a genotype's offsets, with independent draws per fly.
Each fly's record is rotated circularly by its own offset r, uniform on {0, ..., 1439} minutes (never whole days,
which would keep the alignment). The rotated record's minute t holds the original minute (t + r) mod T. The wrap
junction is a break: the original record has no transition into its first minute, so none is counted across the
junction. Each replicate goes through the identical pooling, binning and cosine fit.
What the null keeps: each fly's own rhythm, noise and sample size. What it destroys: phase locking to clock time
shared across flies. The question is "do the flies share a clock-locked modulation of J?"
p = (1 + #{A_null >= A_obs}) / (B + 1).

## Calibration gate (BEFORE the real verdict; synthetic)
200 synthetic datasets matched to per0 DD (its 25 flies, their record lengths and missing-minute masks). Each fly is
simulated from its own constant-rate observed chain (sqw_ref.mc_fit, whole record, one matrix; sim_rows for empty
rows; sqw_ref.simulate with B = I; missing minutes reinserted), so there is no modulation. Dataset d uses
default_rng(21_000 + d) for the simulation and default_rng(22_000 + d) for its 999 rotations. Full test on each.
GATE: the fraction with p < 0.05 must lie in [0.02, 0.09]. If it fails: STOP, record, report back; the test is not
tuned. The real-data verdict is computed only after the gate passes.
Optional, same run: POWER at per0 DD's sample size. 100 datasets of 25 flies with per0 DD's record lengths, simulated
from CS LD's pooled 12-bin chains P_b (the same time-varying chain for every fly, transition into minute t uses its
clock bin), so the modulation has CS LD's observed shape and amplitude. Seeds default_rng(23_000 + d), rotations
default_rng(24_000 + d). Power = fraction with p < 0.05. Stated limitation: no between-fly heterogeneity, so the
power is optimistic.

## Verdict (one run, no re-runs)
- Holm over the three primary genotypes at alpha 0.05; per genotype "clock-locked modulation of J" if Holm p < 0.05.
- per0 DD should be null. If its p < 0.05, the whole test is FLAGGED: every positive must carry the caveat "null
  miscalibrated or residual rhythm in per0 DD", and none may be reported as clean.
- LD positives are light-driven plus clock; this test does not separate the two, and the write-up must say so.
- Narrowest block: each genotype's whole recording, pooled in 2-h clock-time bins. The null is sample-size matched by
  construction (same flies, same lengths).
Figure: the 2-h J profiles with the fitted cosine overlaid and the null's 95% amplitude band.
