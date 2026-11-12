# beam_hmm results ledger (append-only)

Every result, null results included, with date and the commit hash of the code that produced it.

## 2026-11-12 — E0 exploration of week 1 (code 2fa43a6; PRESPEC E0; descriptive, no tests)
Format
- Positions every 0.4 s, not 0.2 s: data lines alternate with all "-1" lines that carry no time stamp
  (1,511,985 data lines). 18 irregular steps (0.6-1.2 s) at t ~ 220-221 ks, 305-307 ks, 488 ks.
- Channels 33-64: constant (4) except monitor-wide flips at t ~ 551.8-552.4 ks -> no flies.
- Channels 17-32 (presumed monitor 2): values 0-16; value 15 essentially never occurs (dead beam or unused code).
- Channels 1-16 (presumed monitor 1): values 13-16 ONLY, in every channel. Encoding unknown; not interpretable
  as position until explained.
Artifacts (monitor 2)
- Clean (>= 95% of changes <= 2 beams per 0.4 s): ch18, 19, 22, 24, 26, 28, 29.
- Flicker/hub artifacts (many single-sample jumps >= 4 beams to one hub value): ch17, 21, 23, 25, 27, 31, 32
  (e.g. ch31: 6508 single-sample spikes, hub value 3).
Silence / death
- ch16 and ch30: <= 4 changes all week (empty or dead from the start). ch3: silent 82.6 h from 73.7 h, then
  brief changes at ~156 h. ch27: silent from 132.8 h. ch32: silent for the first 17.8 h.
- Common event at ~153-157 h (day 7, ~ZT9-13): many channels burst together, then ch6, 7, 19, 22, 25, 29, 31, 32
  go silent to the end of the file (11-14 h). Simultaneous -> handling/transfer, not death. A 12-h silence rule
  applied blindly would call some of these dead; handling can also move dead flies (ch3, ch30).
Descriptive (clean channels 18,19,22,24,26,28,29; hours 0-150)
- Time in runs of no beam change: >= 5 min 71%, 1-5 min 13.5%, 5-60 s 10.6%, < 5 s 4.9%.
- Midday occupancy concentrates at beam 16 in ch18, 26, 29 (partly 19, 28), at low beams in ch22, 24.

## 2026-11-12 — P1 on week 1 (code 3151481; PRESPEC P1, pilot; descriptive)
Inclusion (channels 17-32, week 1)
- Included (8): ch18, 19, 20, 22, 24, 26, 28, 29.
- Artifact-heavy (7; > 5% of changes are >= 4-beam jumps after the spike filter): ch17 (13.9%), 21 (9.3%),
  23 (8.0%), 25 (8.1%), 27 (20.4%), 31 (52.8%), 32 (43.0%).
- Empty/dead at load (1): ch30 (0 changes in the first 24 h).
HANDLING DETECTOR FAILED: the pre-specified rule (>= 50% of included channels above their own 99th
percentile in one 10-min window) found NO window in week 1, although E0 shows a common event at ~153-157 h.
So the whole file after hour 2 was treated as valid, and the 13.6-14.4-h silences that end the file for
ch19, 22, 29 are counted as valid inactivity (they inflate day-7 sleep). Not re-tuned on these data: the
rule will be replaced by the user's transfer times, or by a new dated detector entry tested on weeks 2-9.
Descriptive (included flies; per-fly values, median over minutes; range over 8 flies)
- Sleep per day (whole week, inflated as above): 864-1281 min (median ~1080; 14-21 h). ch20 and ch26 are
  the 1280-min flies and also the least active in the first 24 h (396 and 204 changes).
- Share of active seconds that contain ONLY edge flicker (1-beam A-B-A): S minutes 0.9-10.8% (light),
  1.2-12.9% (dark); Q 0-5.5% / 1.0-8.8%; W 0.9-6.5% / 0.5-7.4%.
- W minutes: median beam changes per minute 3.7-15.5 (light), 3.8-22.7 (dark); longest still run inside a
  W minute median 25.6-60 s (light), 15.6-60 s (dark).
- Per-minute medians are dominated by zeros for S and Q; distribution-level summaries are needed (next).

## 2026-11-12 — D1 gate code smoke test (code a8447bc; toy settings, NOT the gate)
1 synthetic fly per scenario, 4 days, K in {2,3,4}, 1 restart, 60 iterations, 2 bootstrap replicates.
- G1: K* = 4, but the four states were quiescent-END-S, quiescent-BEYOND1-S, quiescent-MIDDLE-S and
  active-LOCAL-W: sleep split by place; the planted BOUNDARY state's best match had Jaccard 0.06.
- G2: K* = 3; a spurious candidate passed C1 and was "promoted".
Reading: with zone in the symbol, persistent position during rest makes place look like hidden state.
Logged, not acted on: the gate (full settings) runs as pre-specified.
