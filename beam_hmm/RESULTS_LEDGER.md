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

## 2026-11-12 — E0 follow-up: jump artifacts have a hub (code 157bf37; week 1, descriptive)
Per channel: share of >= 4-beam jumps going to the single most common target, and share of runs at that
target that return to the pre-jump value. Artifact-heavy: ch17 16 (45%, returns 76%), ch21 5 (44%, 21%),
ch23 3 (44%, 30%), ch25 8 (40%, 45%), ch27 16 (45%, 52%), ch31 3 (47%, 35%), ch32 11 (41%, 55%); median stay
at the hub 1.2-3.8 s. Clean channels: top target 15-39%, returns 0-8%. Reading: an intermittently blocked
beam whose reading overrides the fly's position.

## 2026-11-12 — R1 synthetic check of the P1-B hub repair (code 40b5c3a): FAIL
10 synthetic G1 flies x 7 days with a hub artifact. Hub identified correctly in 10/10 (share 0.50, returns
0.94-0.97); no false hub in the 10 artifact-free flies (0 samples changed).
S/Q/W minute agreement with the artifact-free fly: before repair 0.21-0.33, after 0.61-0.84 (need >= 0.95).
Window-symbol agreement: before 0.86-0.88, after 0.93-0.97 (need >= 0.95).
Cause: the rule only repairs hub runs entered by a >= 4-beam jump; when the fly rests within 3 beams of the
hub, blocking episodes look like small steps, are kept, and each one breaks a sleep bout.
Per P1-B, repaired data are NOT used; the rule is changed in a new entry (P1-C).

## 2026-11-12 — R1 re-run under P1-C (code 17d5de1): PASS
Hub found in 10/10 synthetic flies; S/Q/W agreement after repair 0.969-0.993 (before 0.21-0.33); window
symbols 0.991-0.995 (before 0.86-0.88); artifact-free flies changed: 0 samples in 10/10.

## 2026-11-12 — P1-C hub repair applied to week 1 (code 17d5de1; QC only)
Dead-hub rule fires in exactly the 7 artifact-labelled flies and in none of the 8 others (their top-target
shares 0.15-0.39 with returns 0.00-0.25 do not meet both conditions).
Samples replaced / >= 4-beam jumps remaining (share of changes): ch17 8.0% / 3.1%; ch21 12.1% / 2.4%;
ch23 21.8% / 3.1%; ch25 11.7% / 2.5%; ch31 29.1% / 11.3%; ch32 9.9% / 18.1%; ch27 97.6% / 27.5%.
Clean flies for comparison: 0.1-3.7% jumps >= 4.
ch27 is essentially unreadable (hub 16 holds 97.6% of its readings). ch31 and ch32 keep a second artifact
source (e.g. ch32's second target, 16, takes 39% of its big jumps).

## 2026-11-13 — D1-B synthetic gate (per-fly results from code 752e34f; REPORT.txt written at 9790253, gate_d1.py
## and d1.py identical between them; run on the user's machine, pushed as 17a458c): FAIL
G1 (planted BOUNDARY state): K* = [5,6,6,6,7,5,4,5,5,4,6,4,6,4,4,4,5,4,5,5] -> K* = 4 in 7/20 (need 14). Planted state:
candidate 20/20, C1 19/20, C2 20/20, C3 17/20, C4 0/20 (expected: flicker-defined), promoted 17/20 (need 16);
matched Jaccard 0.66-0.75, so the secondary (D1-A2) count is also 17/20.
G2 (no extra state): K* = 3 in 11/20, 4 in 6/20, 5 in 3/20; spurious C1 18/20, promoted 18/20 (allowed 2 and 1).
Post hoc reading of the per-fly fits (state types and dwells in the JSONs; not a pre-specified analysis):
- Extra states: K* > true K in 13/20 G1 flies; the surplus states have expected dwell 10-14 s (one window) in
  most of them, e.g. ACTIVE/LOCAL/W 10 s, ACTIVE/TRAVEL/W 10 s, and add 0.0001-0.0006 nats/window of held-out
  logL (a real state adds ~0.017). Consistent with the windowing explanation behind D1-D.
- The G2 rule itself is unpassable: with K* = 3 the three fitted states ARE the true ones (QUIESCENT/S ~1700-
  2400 s, ACTIVE/LOCAL/W ~120-160 s, ACTIVE/TRAVEL/W ~40-60 s), but true QUIET and TRAVEL are both scored W,
  so they "split a label", become candidates, pass C1 against K = 2, and count as spurious. A gate design
  error (mine), fixed for D1-D in entry D1-D2. G2 also fails on genuine over-splitting in its 9 flies with
  K* > 3, so D1-B fails under either reading.
- C1's bootstrap thresholds were ~0 (null simulated exactly from the window-level HMM), so tiny gains passed.
- The T1 gate's TG3 uses the same "any promoted candidate" rule and inherits the same flaw; TG1/TG2 do not.

## 2026-11-13 — T1 gate INTERIM sanity check (30/40 flies, pushed as 10b7f8d; code ac6007d) — not the verdict
Looked at to check output sanity before the gate finishes; the pass rules are fixed and unchanged.
- Sleep-state matching works: Jaccard with true SLEEP 0.98-0.99 in all 18 TG1 flies.
- TG1: a clock model (M1-M3) chosen over M0 in 13/18 (needs >= 16/20 at the end; with 2 flies left it cannot
  reach 16, so TG1 will FAIL on selection). The 5 flies choosing M0 had held-out gains of M1-M3 of
  -0.2 to -11e-4 nats/window; 3 of them had K* = 4 (an extra, over-split state). Where a clock model was chosen,
  the fitted/true peak-to-trough ratio was 0.80-1.64.
- TG1 peak criterion is UNINFORMATIVE as specified (design weakness, mine): the true total sleep-exit
  profile = light/dark base x modulation peaks at ZT 12.0 (the last light window, where the higher light-phase
  base meets the rising modulation), so any fit with the step and a roughly right shape scores 0 h, and any
  M0 fit scores 12 h. It re-measures "clock model chosen", not peak placement. Not changed (gate running).
- TG2: M0 chosen in 11/12 (allowed <= 3/20 false selections); the single false one (fly 15, M1) had a gain of
  +0.05e-4 nats/window.
