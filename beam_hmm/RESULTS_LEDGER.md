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
