# Design note for the repeat experiment (2026-11-15)

Planning document, not a result. Numbers come from the pilot (one run, CS males, LD 12:12, channels 17-32) via
RESULTS_LEDGER.md and `design_power.py` (commit 6633f48; output `results/design/power.json`).
Everything here should be turned into PRESPEC entries **before** the new data are seen.

## 1. How many flies per question

Method: the pre-specified test for each A-measure is a one-sample Wilcoxon signed-rank test across flies, with Holm
correction. I simulated it at the Holm worst case (alpha/m), with per-fly values drawn from Normal(d, 1).
d = pilot median / pilot fly-to-fly SD. Each measure is computed at two effect sizes:
- the pilot effect;
- half the pilot effect. Effects estimated from about 10 flies, and picked because they were significant, are
  usually inflated (Button et al. 2013, *Nat Rev Neurosci* 14:365). Plan on the half-effect column.

| Question (pilot status) | d (pilot) | Flies for 80% power at pilot d | Flies for 80% power at d/2 |
|---|---|---|---|
| A1 sleep at tube ends by day (SUPPORTED) | 1.99 | 8 | 15 |
| A4 slower walking at night (SUPPORTED) | 1.16 | 12 | 40 |
| A5 M1 sleep rises with age (nominal) | 1.13 | 15 | 40 |
| A5 M3 end-resting rises with age (nominal) | 0.51 | 50 | >100 |
| A5 M4 walking speed rises with age (nominal) | 0.36 | 100 | >100 |
| A5 M2 activity per waking hour (n.s.) | 0.27 | >100 | >100 |
| A2 pauses before sleep, illustrative effect 0.03 (n.s.) | 0.43 | 60 | >100 |
| T1-C1 cohort clock modulation (SUPPORTED, p 0.001 at 11 flies; LR 1012 vs null median 81) | n/a | already saturated at 11 | — |

Each entry is the number of **analysable** flies. "Analysable" means a fly in the primary set for that test.

**Yield in the pilot** (from 16 readable channels):
- 12 flies analysable in week 1;
- 11 in weeks 2-3;
- 9 with 3 or more qualifying weeks for the ageing slopes (weeks 2-5).

So about 70% of loaded flies are usable for single-week tests, and about 55% for ageing. LoadingCheck should
recover some of the losses: the 14 hidden fly-weeks were all caught in the first 2 h. But that gain has not been
measured, so plan on the pilot yields.

**Recommendation:**
- **Minimum: about 30 loaded flies.** Gives about 20 analysable flies. That covers A1 and T1-C1 even at half the
  pilot effect, but covers A4 and M1 only if the pilot effects are real at full size.
- **Recommended: about 70 loaded flies.** Gives about 50 analysable for single-week tests and about 40 with 3 or
  more weeks. That covers A1, A4 and sleep-versus-age (M1) at half the pilot effect. M3 and M4 are out of reach at any realistic n unless slopes become more precise per
  fly. Each slope's SD shrinks with more weeks per fly and more valid days per week. Fewer transfers would help,
  as would a shorter first-day exclusion if it can be justified *in advance*.
- **A2 / A3 (small or absent effects):** treat as exploratory. Detecting effects this small needs 60-100+ flies.

**Caveats:**
- d uses the median, not the mean.
- Values are assumed normal.
- Holm worst case is assumed for every test, which is conservative.
- The slope SDs depend on how many weeks each fly contributes. The pilot's flies contributed 3-4 weeks.

## 2. Data-quality protocol (each item answers a pilot failure)

1. **LoadingCheck within 2 h of every loading and transfer**, and re-seat any tube it flags. In the pilot it caught
   all 14 fly-weeks in which a live fly was hidden behind a blocked beam 16 (about 1 fly-week in 10). It also caught
   about half of the blocked beams directly, and about 70% with any flag. Freeze its thresholds now; changing them
   needs a new PRESPEC entry.
2. **Record which end the food is at, for every tube.** This tests whether beam blockages come from the food or the
   plug. It also lets A1 ("flies nap at the tube ends") be split into food end vs plug end, which is your hypothesis
   that flies nap near food. The pilot cannot test this.
3. **Record clock times:** lights-on and lights-off by the incubator, the recording computer's clock offset at each
   transfer, and the exact loading and transfer times. In the pilot, the phase had to be inferred from the flies'
   own lights-off rise (P1-E), because files started 30-40 min off and week 6 needed a hand check.
4. **Check every channel before the run**, with a 1-h empty or test recording. In the pilot only 16 of 64 channels
   were readable.
5. **Keep the narrow-range rule** (fewer than 8 beams spanned in a fly-week means exclude), and keep the P1-C
   dead-hub repair and the P1-F death rule unchanged. They are already pre-specified, so they carry over.

## 3. Design elements for the open questions

- **A DD segment** (e.g. 4-5 days of DD after week 2, before ageing confounds it). In LD, T1-C1 cannot separate a
  clock from lighting effects. It also cannot settle the A2b direction question: the pilot (CS, LD) shows a slow
  wake-up, while the earlier per0 / CS DD data show a wind-down into sleep. A CS DD segment in the same monitors
  tests whether the condition is the cause.
- **Optional genotype contrast:** a per0 group run alongside the CS flies, with matched ages and the same monitors,
  would separate genotype from cohort. It is optional because it doubles the fly count.
- **Approach to death (Aim 3):** in the pilot, few valid days fell close to death. Each 7-day file loses its first
  24 h, and transfers add gaps. The options are fewer transfers (if food lasts), or a pre-specified shorter
  exclusion after transfer, for example 2-6 h, justified from the transfer-response time course in the existing
  pilot data. Fix that choice before the new data arrive.
- **Food end and siesta position:** with item 2 recorded, A1 can become a three-way comparison: food end vs plug
  end vs middle.

## 4. Compute budget for the model work (Track B)

The B1 synthetic gate takes about 15-18 core-hours per 14-day fly on the user's PC. Most of that goes on the
19-dataset bootstrap behind each fly's state-count decision.
- **If B1 passes:** fitting about 40-70 real flies × several weeks with the same procedure would take thousands of
  core-hours. That needs a cluster, or a cheaper pre-specified procedure such as fixed K plus a stability check,
  or fewer bootstrap datasets with a sequential stopping rule.
- **If B1 fails:** model work moves to fixed-K and semi-Markov descriptions, which are much cheaper.

Either way, decide this before the new data arrive.

## 5. What to pre-specify before the data arrive (checklist)

- [ ] Test weeks and primary sets for A1, A4 and A5 M1 (replication of the pilot tests, same code).
- [ ] T1-C1 on the new LD weeks, plus the DD segment as its own pre-specified test.
- [ ] A2b on the DD segment (direction question).
- [ ] Food-end split of A1.
- [ ] LoadingCheck thresholds (frozen) and the re-seat rule.
- [ ] Transfer-exclusion window.
- [ ] Track B procedure, chosen after the B1 verdict.
