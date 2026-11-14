# beam_hmm plan (2026-11-14, user decision: option 3)

Status of the data: a PILOT (channels 17-32 of one 9-week CS run; 12 usable flies in week 1). A larger,
repeated experiment is expected in the coming months. Goal of the pilot phase (user, 2026-11-14): explore
broadly enough that, once there is enough data, we can settle quickly on the best approach. So every
track below ends in a DECISION RULE for the large dataset, not in a pilot claim.

## Where the gates stand (all synthetic; details in RESULTS_LEDGER.md)
| gate | what it tested | verdict | lesson |
|---|---|---|---|
| D1-B | state discovery, 10-s windows | FAIL | windows create transitional states (dwell ~10 s) |
| D1-D | state discovery, 0.4-s steps | FAIL | degeneracy rule wrong at this level; serial correlation within states makes spurious states as large as a real one |
| T1 | per-fly clock selection | FAIL (TG2 pass) | does not invent rhythms; too little power fly by fly |
| T1-B | cohort clock test on sleep entry/exit | queued | -- |

## Track A — targeted descriptive analysis of the pilot (now)
Model-free (or model-light) measures, each with flies as the unit and a sample-size-matched null where a test
is made. Candidates named so far (E0, P1, chat 2026-11-12):
1. Where flies rest: share of sleep spent at a tube end vs the middle, and beyond beam 1, by ZT; compared
   with where the same fly is awake (within-fly contrast, so the unknown food end does not matter).
2. Pauses that never reach the 5-min criterion (5 s - 5 min) inside wake: rate, durations, place, ZT.
3. Edge flicker: how much "activity" is one-beam flicker; how sleep amount and bout structure change if
   flicker is not counted as movement (a sensitivity analysis of the sleep criterion itself).
4. Walking: speed and run length (beams per 0.4 s; distance before a stop), by state and ZT.
5. Ageing (when weeks 2-9 are in): each of 1-4 per fly per week, aligned by age AND by days before death.
Each becomes a PRESPEC (A1, A2, ...) before it runs. What A gives the large dataset: the measures, their
fly-to-fly spread (for power calculations), and which features are worth modelling.

## Track B — models, gated on synthetic data (in parallel, on the user's machine)
- B1 (= D1-E): direction-aware autoregressive HMM at 0.4 s. Each step coded relative to the previous one
  (stay / same direction / reversal / >= 2 beams), emissions conditional on the previous step, so travel
  persistence and flicker pairs no longer force extra states. Degeneracy defined on emissions AND dynamics.
  Same G1/G2 truths and D1-D2 scoring. NEW for breadth: the gate is run at 7 and 14 days per fly, so it
  also says how much data per fly the method needs.
- T1-B: cohort clock test (queued; K fixed at the truth, so independent of B1).
What B gives the large dataset: a validated model family (or the knowledge that discovery by likelihood
does not work at this resolution) and the data needed per fly.

## Decision table for the large dataset (to be filled as gates and analyses report)
| question | candidate approaches | chosen when |
|---|---|---|
| sub-minute states exist? | B1 discovery; A2/A3 targeted measures | B1 gate passes -> B1 on the new data; else A-measures only |
| clock control of sleep entry/exit | T1-B cohort test | T1-B gate passes |
| ageing | A-measures per week; B1/T1-B per week | after A5 shows which measures carry fly-to-fly signal |
| how many flies / days needed | gate runs at 7 vs 14 days; A-measure spread | before the repeat experiment is designed |

## Needed from the user
- Weeks 2-9 (channels 17-32 are enough).
- Transfer times (or a go-ahead for a new handling detector, tested on weeks 2-9).

## Open question: direction of the seconds-scale wake-bout asymmetry (A2b), recorded 2026-11-14
Pilot (CS LD, this run): first ~20-25 s after waking LESS active than the last ~25 s before sleep (secondary setting).
Earlier project (per0 DD, CS DD; and CS LD descriptively): the reverse (wind-down into sleep, abrupt waking).
Ruled out (ledger): porting error (exact reproduction of per0 DD), detector type (both multibeam), activity conversion
(lab conversion gives the same pilot direction), sleep scoring (identical second for second).
Agreed with the user: the difference is in the flies or the conditions -- cohort, age, genotype, LD vs DD -- or in the
one remaining method difference (the pilot excludes bouts within 1 h of a light switch; the earlier CS LD block did
not). Earlier CS LD already showed the earlier direction, so LD alone is an unlikely explanation.
For the repeat experiment: CS in LD AND DD in the same cohort, ages recorded, A2b's 5-s setting pre-specified as
primary, and the light-switch exclusion reported both ways.
