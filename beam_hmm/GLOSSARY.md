# beam_hmm glossary

Plain-language definitions of the terms used in the beam_hmm tools, printouts and results. Each entry gives
what the term means physically, the exact rule, where the rule is fixed (PRESPEC entry), and which function
reports it. When this file and a PRESPEC entry disagree, the PRESPEC entry is the definition; please report
the discrepancy.

Functions: MATLAB in `beam_hmm/matlab/` (`help LoadingCheck`, `help RawToActivity1s`, `help BeamFingerprint`);
Python twins in `labconv.py`, `pipeline.py`, `p1e.py`.

---

## 1. The recording

**Channel / tube / fly No.** One tube in the monitor, one fly per tube. The fly number on the survival sheet
equals the channel number. In the 2026 pilot only channels 17-32 were readable.

**Beam, position.** Each tube is crossed by 16 infrared beams. The position is the beam that is currently
broken: 1 to 16. Beams 1, 2, 15 and 16 are the **end beams**.

**Position 0 ("no beam broken").** No beam is interrupted. In the pilot, 0 was always entered from beam 1, so a
live fly reading 0 is beyond beam 1, in the detection-limited end of the tube. **In a correctly oriented tube
beam 1 is the beam closest to the food, so 0 means the fly is at the food and very likely eating** (user,
2026-11-16). It is a proxy for "at the food", not a measure of ingestion; in a tube mounted the other way
round, 0 is the plug end. Record the orientation of every tube. Dead flies often read 0 as well.

**AtZero / EntriesToZero (feeding proxy).** `RawToActivity1s` returns `out.AtZero`: for every second, the share
of its readings at 0 (0..1), after the same repairs as the activity; optionally written to `AtZeroFolder` (never
the activity folder, which `process_LocomToSleep` reads in full). Its QC table gives `ShareAtZero` and
`EntriesToZeroPerDay`. `LoadingCheck` gives `EntriesToZero` in the window and lists active tubes with no visit
to 0 (information, not a flag: possible reversed tube, or food or plug filling the space beyond beam 1). None of
these changes the activity counts or the sleep scoring; they let a sleep bout spent at 0 be flagged. `RawToActivity1s` keeps 0 as a position by default (`'KeepBeyondBeam1', true`), so crossing beam 1
counts as movement. The lab's `PreProcess.m` folds 0 into the previous reading instead.

**Beam 15 (dead).** In this monitor beam 15 never reports. A fly walking from 14 to 16 therefore shows up as
a 2-beam step. This is a hardware fact, not something the pipeline repairs.

**Raw file format.** 3 header lines; then per line 64 positions and a time column. Every other line is
64 x -1 with no time; these lines are dropped (`'beam_hmm'` mode). Source: `PROTOCOL.md`.

**Sampling step (0.4 s).** Measured from the time column after dropping the -1 lines. The lab's older code
treated every line, -1 lines included, as 0.2 s. `RawToActivity1s` measures the step and reports it.

**Position change.** Two consecutive timed readings that differ. The basic count behind most flags.

**Big jump.** A change of 4 beams or more between consecutive readings (0.4 s apart). A walking fly moves
about 1-2 beams/s on average (A4 medians), so a big jump in one step is physically implausible and almost
always means a detector artifact.

---

## 2. Artifacts and their repair (applied before any analysis)

**Spike.** A single reading B that jumps 4 or more beams away and comes straight back: A-B-A with
|B - A| >= 4. **Spike filter:** B is set to A, in one pass. Defined in PRESPEC P1; applied by
`RawToActivity1s` ('beam_hmm' mode) and `pipeline.spike_filter`.

**Hub (blocked beam).** A beam whose reading intermittently overrides the fly's position, e.g. because food,
the cotton plug, a droplet or debris partly covers it. The trace then jumps to that beam and back,
repeatedly. A hub is identified per channel and recording (or per window in LoadingCheck) by three numbers:
- **Hub value**: the beam that is the most frequent target of big jumps.
- **Hub share**: the fraction of all big jumps that land on the hub value. Must be >= 30%.
- **Hub returns**: of the visits to the hub that begin with a big jump, the fraction that end by returning to
  the reading before the jump. Must be >= 20%. (A fly that really walks to a beam does not keep returning to
  exactly where it was; a beam flickering on and off over a stationary fly does.)

Defined in PRESPEC P1-B (identification). A fly that rests on one beam for a long time is NOT a hub: only the
jump signature counts.

**Dead-hub repair ("hub declared dead").** If a hub qualifies, that beam is treated as dead, like beam 15:
every reading of it is replaced by the last reading that was not the hub. PRESPEC P1-C (replaced P1-B's first
repair rule, which failed its synthetic check R1). Cost, stated in P1-C: time the fly genuinely spends at the
hub beam, and movement onto it, become invisible. Validated on synthetic flies with a planted blocked beam:
sleep/quiet/active minutes agree 97-99% with the artifact-free truth (P1-C, R1).
Reported by `RawToActivity1s` as `HUB k declared dead (x% of readings replaced)`; in 'legacy' mode as
`HUB k DETECTED, NOT repaired` (its flicker then counts as movement).

**End-beam hub.** A hub at beam 1, 2, 15 or 16. Points to food or the plug at that end. The fly's position at
that end of the tube cannot be known. Flag: `END-BEAM HUB` (`RawToActivity1s`).

**Narrow range.** The positions span fewer than 8 beams (max - min < 8) although the fly moves. Part of the
tube is obstructed, or the channel reads only a few beams. Exclusion rule fixed in A2-A4 (2026-11-14) and in
P1-F (as UNREADABLE). LoadingCheck flags it as `NARROW RANGE` when there are >= 50 changes in the window.

**Artifact-heavy fly-week (P1-D).** More than 5% of the position changes are big jumps AFTER spike filter and
hub repair. Excluded from pooled analyses and reported separately.

**Flicker (A3).** A single-sample 1-beam out-and-back (A-B-A with |B - A| = 1). Not removed by the spike
filter (that needs >= 4 beams). In the pilot its rate did not differ by day and night, which points to a
beam-boundary effect (a fly sitting on the edge of a beam) rather than behaviour; the 5-min sleep measure
changes by only 1-5 min/day when flickers are ignored (ledger, A3).

---

## 3. LoadingCheck flags (first 2 h after loading or re-seating)

All computed on the window only (default the first 2 h; `'StartHour'` moves it). Thresholds are fixed counts
and were not tuned on the pilot. Source: `LoadingCheck.m` header; Python twin `labconv.loading_check`.

| Flag | Rule in the window | Physical meaning | Action |
|---|---|---|---|
| `BLOCKED BEAM k (end / mid-tube)` | >= 10 big jumps, hub share >= 30%, hub returns >= 20% | beam k blocked SOMETIMES; the fly is visible most of the time | data usable after dead-hub repair; re-seat if easy |
| `FROZEN AT END BEAM k` | one value in >= 99% of readings, at beam 1, 2, 15 or 16 | blocked ALL the time (or a dead fly, or a fly resting at the end beam): a live fly can be completely hidden | re-seat now; check plug and food |
| `FROZEN AT BEAM k` | the same at a mid-tube beam | stuck reading, dead fly, or long rest | check the tube and fly |
| `FROZEN AT 0` | the same at 0 | fly beyond beam 1 the whole window: at the food (correct orientation; feeding or resting there), or dead | check the fly |
| `NARROW RANGE` | >= 50 changes but < 8 beams spanned | part of the tube obstructed | check how far food and plug reach in |
| `LOW ACTIVITY` | < 20 position changes | dead or stuck fly, blocked tube, or a sleeping fly | check the fly |

**BLOCKED vs FROZEN.** BLOCKED = blocked sometimes, data recoverable. FROZEN = blocked all the time, data
lost while it lasts. Both on one channel: treat as FROZEN.

**Re-check labels (`'Previous'`).** Compared with an earlier check of the same tubes:
RESOLVED (flagged before, OK now), STILL FLAGGED, NEW FLAG (OK before, flagged now -- e.g. a disturbed
neighbour), OK (OK before and now), NOT IN PREVIOUS (channel not checked before). RESOLVED means "not seen
in this window" only.

**Limits (pilot, 144 fly-weeks).** The 2-h screen flagged all 14 live flies hidden behind a blocked beam 16,
and ~half of the intermittent blockages found later in the full week (~70% with any flag); 1 h gave the same
catches. A window during the siesta or at night flags sleeping flies as LOW ACTIVITY or FROZEN: re-check soon
after handling. Ledger 2026-11-14 and 2026-11-15/16.

---

## 4. Readability, outages and death (whole fly-weeks)

**EMPTY/DEAD (`RawToActivity1s` status).** Fewer than 100 position changes in the first 24 h of the file.

**DEAD-LIKE / UNREADABLE / READABLE (P1-F).** Per fly-week, after spike filter and hub repair:
DEAD-LIKE = fewer than 100 position changes in the whole fly-week; UNREADABLE = >= 100 changes but the hub
repair replaced more than 50% of the samples, or the positions span < 8 beams; READABLE otherwise.

**Outage.** A DEAD-LIKE (or frozen) fly-week followed by an active READABLE week. Not a death. In the pilot,
14 fly-weeks read a constant 16 for the whole week while the fly was alive (hidden behind a blocked beam 16).

**Death (P1-F).** A fly is dead from the start of the last >= 12-h silence in its last READABLE fly-week, if
every later fly-week is DEAD-LIKE. If a later fly-week is UNREADABLE, the death is UNDETERMINED and the fly is
excluded from that point (never counted as dead). The 24 h before a death are excluded (the lab's rule).

---

## 5. Time

**ZT (zeitgeber time).** Hours since lights-on; ZT 0-12 light, ZT 12-24 dark (LD 12:12).

**Behavioural phase stamp (P1-E).** Because the recording computer's clock drifted from the incubator's,
ZT is taken from the flies: pooled activity in 10-min bins, folded at 24 h; LIGHTS-OFF = the bin with the
largest rise; LIGHTS-ON = 12 h earlier. Guard: an estimate more than 1.5 h from the file start is flagged and
checked by hand (week 6 of the pilot was). MATLAB counterpart: `EstimatePhaseStamp` (in the lab's folder).

**First-cycle exclusion (P1-E).** The first 24 h of every file are excluded (acclimation after loading).

**Transfer.** The weekly move to a fresh-food tube; file boundaries are transfers.

---

## 6. Sleep and wake scoring

**Activity per second.** Beams travelled in each second (sum of |step| over the position changes whose time
falls in that second). `RawToActivity1s` output; the lab's `multibeam_locomotion` unit.

**Sleep bout.** >= 5 min (300 s) without movement (activity 0 in every second), with the lab's bridging: an
interior wake run of <= 1 s inside sleep is turned into sleep (`ironout`). Identical to the lab's
`sleeptrace` + `ironout`, second for second (ledger: 0 of 28M s differ).

**Wake bout.** The time between two sleep bouts.

**S / Q / W minutes (P1).** Per minute: S (sleep) if >= 30 s of the minute lie in a sleep bout; otherwise
Q (quiet wake) if the share of active seconds is below half of the fly's reference activity for that phase
(the median activity of the middle 5 min of its wake bouts >= 10 min); otherwise W (active wake).

**Pause.** A run without any position change lasting 30 s to < 5 min, outside sleep bouts (A1, A2).

**Sleep onset / wake-up.** The first second of a sleep bout / the first second after it.

---

## 7. Measures used in the Track A analyses

**END / MIDDLE / BEYOND1 zones.** END = position 0, 1, 2, 14, 15 or 16; BEYOND1 = position 0; MIDDLE = the rest
(D1, A1).

**D (A1).** Share of sleep samples at END minus share of pause samples at END, per fly and light/dark block.
D > 0: the fly rests at the tube ends more when asleep than in short pauses.

**Walking bout, speed (A4).** A maximal run of position changes no more than 2 s apart with path length >= 3
beams. Speed = path length / (time from first to last change + 0.4 s), in beams/s.

**Pause fraction, PRE / MID / POST (A2).** Share of a 10-min window inside pauses; PRE = the 10 min before a
sleep onset, MID = the middle 10 min of the wake bout, POST = the first 10 min after waking.

---

## 8. Model terms (Track B)

**HMM.** Hidden Markov model: the fly is assumed to switch among a few hidden states, each with its own
statistics of what is observed.

**B1 symbols.** Each 0.4-s step coded relative to the last move: STAY, SAME (1 beam, same direction), REV
(1 beam, reversed), BIG (>= 2 beams). B1's states predict the next symbol from the previous one
(autoregressive emissions).

**K\*.** The number of hidden states chosen for a fly by held-out likelihood (fit on half the days, score on
the other half).

**Degenerate states (B1).** Two states with nearly the same emissions (< 0.05 bits apart) AND similar dwell
times (within a factor 1.5). Fits with degenerate states are not eligible for K\*.

**C1-C4 (promotion criteria).** C1: held-out gain beats a simulated null from the K\*-1 model (95th
percentile of 19 replicates); C2: mean dwell >= 30 s; C3: the state is found in both halves of the data;
C4: survives flicker removal and is not mainly the fly beyond beam 1. A state that passes is "promoted".

**Spurious state (D1-D2 scoring).** A fitted state whose best overlap (Jaccard, over samples) with every true
simulated state is < 0.5.

**Gate.** A synthetic test with known truth that a method must pass before it is used on real flies; pass
rules fixed in advance (G1 positive control with a planted state, G2 negative control without one).

---

## 9. Process words

**PRESPEC.** `PRESPEC.md`: the analysis plan, committed before the analysis runs, so its timing is checkable.

**Ledger.** `RESULTS_LEDGER.md`: every result, null results included, with the date and commit hash;
append-only (corrections are appended, never edited in).

**Post hoc / descriptive.** Computed after seeing the data, or not part of a pre-specified test; reported but
not a finding.

**Primary / sensitivity set.** Primary = fly-weeks passing the inclusion rules (P1-D); sensitivity = flies never
artifact-heavy, reported beside it.
