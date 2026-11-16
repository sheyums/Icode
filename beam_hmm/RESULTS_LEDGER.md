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

## 2026-11-13 — T1 synthetic gate, TG1 + TG2 (code ac6007d; run on the user's machine, pushed as 2622927): FAIL
TG3 deferred (T1-A). TG2_fly09 did not finish (~13.4 h, skipped at the user's request; note in
results/t1_gate/NOTE_TG2_fly09_skipped.txt); neither verdict depends on it.
TG1 (planted morning/evening sleep-exit rhythm, peak-to-trough 3): clock model chosen in 13/20 (need >= 16) -> FAIL.
  M* = [2,2,2,0,2,1,0,0,1,2,1,2,1,2,0,0,0,0,1,2]. Peak within 2 h: 13/20 (the same 13: the criterion is
  uninformative as specified, see the interim entry). Median fitted/true peak-to-trough 0.95 -> within [0.67, 1.5].
TG2 (no rhythm beyond the step): clock model chosen in 1/19 (allowed <= 3) -> PASS.
Reading: the clock models do not invent rhythms and recover the rhythm's size when chosen, but per-fly
held-out selection detects a 3-fold sleep-exit rhythm in only ~65% of flies at 14 days.
POST HOC (computed after the verdict, from the per-fly held-out gains; not pre-specified): resampling 12-fly
cohorts from these flies, T1's own across-fly decision rule (one-sided sign test on best-of-M1/M2 minus M0)
rejects in 15% of cohorts with the rhythm present (Wilcoxon signed-rank 22%) and 0% without it. So at the
pilot's size, the pre-specified real-data test has low power for a rhythm of this size.

## 2026-11-14 — D1-D synthetic gate (code 2f10f70, D1-D2 scoring; run on the user's machine, pushed as 81ff31a): FAIL
G1: K* = 2 in 20/20 (need K* = 4 in >= 14/20); planted state promoted 0/20 (need >= 16) -> FAIL.
G2: K* = 2 in 20/20; spurious C1 0/20, promoted 0/20 -> PASS, but only trivially: a 2-state fit cannot hold
a spurious state (the user's observation), so G2's pass is uninformative.
Diagnosis (post hoc, from the per-fly JSONs):
1 The degeneracy rule carried over from D1 (two emission rows with JS < 0.05 bits) flagged EVERY fit with
  K >= 3 in all 40 flies, so K* = 2 was forced. At the sample level the true states differ mostly in their
  dynamics, not their emissions: SLEEP (1,0,0,0) and QUIET (0.97,0.03,0,0) are only ~0.015 bits apart. The
  0.05-bit threshold was fixed for 10-s window symbols and is wrong at this level (my error in D1-D).
2 Ignoring the degeneracy flag, held-out logL would pick K = 7-8 (the cap) in 38/40 flies. Gain of adding
  the k-th state, median over flies (1e-4 nats/sample):
     K=3: G1 52.9, G2 59.5 | K=4: G1 2.73, G2 2.77 | K=5: G1 2.63, G2 1.14 | K=6-8: G1 2.0-2.5, G2 0.04-0.76
  The REAL 4th state in G1 (the planted boundary rest) adds 2.73; a SPURIOUS 4th state in G2 adds 2.77.
  So D1-D item 8's pre-stated fallback (a minimum held-out gain per added state) CANNOT separate them at
  this resolution: the misspecification of the step-symbol HMM (TRAVEL's direction persistence and the
  out-and-back pairing of a flicker are serial correlations within a state, which a conditionally
  independent emission model can only absorb with extra states) produces gains of the same size as the
  planted state.

## 2026-11-14 — A1 on week 1 (code 081479e; DESCRIPTIVE ONLY, week 1 was viewed before A1 was written)
Window: hours 2-150 of week 1 (first 2 h dropped; stops before the hour-153 handling event seen in E0; chosen
ad hoc for this descriptive run because transfer times are not yet known). Primary set: 12 flies.
D = (share of SLEEP stillness at a tube end) - (share of PAUSE stillness, 30 s - 5 min, at a tube end):
- light: D > 0 in 10/12 flies; median +0.32, mean +0.26, SD 0.19; day-bootstrap 95% intervals exclude 0
  (above) in 9 flies, below 0 in none.
- dark: D > 0 in 6/12; median -0.01, mean -0.04, SD 0.16; intervals above 0 in 1 fly (ch28), below 0 in 5.
Reading (post hoc, NOT a test): by day, flies sleep at the tube ends more than they pause there; at night not,
and several flies sleep away from the ends. The pre-specified TEST is on weeks 2-3 (A1).
For planning: with SD ~0.19 and a mean difference ~0.26, the light-phase contrast would need few flies; the
dark-phase one has no consistent sign.
    CORRECTION (same day): the dark-phase summary above was typed before the summary script ran. From
    results/a1_week1.json: dark D > 0 in 5/12 flies (not 6), median -0.028 (not -0.01), mean -0.052 (not
    -0.04), SD 0.149 (not 0.16). Light values confirmed: 10/12, median +0.321, mean +0.258, SD 0.186.

## 2026-11-14 — weeks 2-5 format and timing QC (code: explore/load_week.py at 6628d30 + inline check)
All four files: same format as week 1 (3 header lines, 64 channels + time, alternating all -1 lines), step
0.4 s with a few irregular steps (largest: one 60.6-s gap in week 3, one 13.6-s gap in week 4).
Phase check (pooled 10-min activity, days 2-6, largest rise and largest rise >= 6 h away): see PROTOCOL.md;
lights apparently on ~05:06 and off ~17:00 by the recording computer's clock in every week; files 1, 2, 5
start ~30-40 min before lights-on. No A1 measure was computed on weeks 2-3.

## 2026-11-14 — A1 PRIMARY TEST, weeks 2-3 (code a427580; PRESPEC A1 + P1-E; weeks 4-5 used only for the death rule)
Phase stamps (P1-E rule 1): lights-on 0.67 h (wk2), 0.00 h (wk3), 23.83 h (wk4), 0.33 h (wk5) into the file; none
flagged. Deaths (P1-E rule 4): ch27 dead from the start of week 2; ch30 from week 4; ch32 possibly dead in week 5.
PRIMARY set, 12 flies (valid days 5.4-10.9 per fly):
- LIGHT: D > 0 in 11/12, median D +0.293, SD 0.147; Wilcoxon two-sided p 0.0033 (normal approximation, because
  one D = 0 exactly), Holm p 0.0067 -> SUPPORTED: by day, flies rest at the tube ends more when asleep than in
  short pauses awake.
- DARK: D > 0 in 4/12, median -0.029, SD 0.136; p 1.0, Holm 1.0 -> not supported, and not the opposite.
SENSITIVITY set: 2 flies (more channels exceed the raw-jump threshold in weeks 2-3 than in week 1): untestable.
POST HOC notes (checked after the result; not part of the decision):
- ch25 (the D = 0 fly) reads only beams 14 and 16 in week 2 (a narrow-range channel, like monitor 1 in week 1),
  so every sample is "at an end" and D = 0 by construction. P1-D has no narrow-range rule (BeamFingerprint
  flags it). Excluding ch25: 11/11 positive, exact two-sided p 0.00098. The conclusion does not depend on it.
- Week 1 (descriptive, before A1 was written) showed the same pattern: light 10/12 positive, median +0.32.
- For the repeat experiment: fly-to-fly SD of D ~0.15 in the light phase; a NARROW-RANGE exclusion rule should
  be pre-specified (range < 8 beams in a fly-week).

## 2026-11-14 — weeks 6-9 format, timing and phase QC (code a427580 p1e.py + load_week.py)
See PROTOCOL.md "Weeks 6-9". Week 6 fails the P1-E phase guard (12-h ambiguity); weeks 7-9 pass (week 9 at
1.50 h, the limit, with 3 flies). Live channels 12/9/4/3 in weeks 6-9.

## 2026-11-14 — A2-A4 PRIMARY TESTS, weeks 2-3 (code 3fb5405; PRESPEC A2-A4 + P1-E; weeks 1-9 supplied for the death rule)
Primary set 11 flies (ch17-24, 26, 28, 29; ch25 excluded by the new narrow-range rule, as intended).
Sensitivity set 2 flies: not tested.
A2 pauses around sleep (Holm over 2):
- PRE - MID pause fraction: median +0.001, 6/11 positive, SD 0.070, exact p 0.966, Holm 1.0 -> NOT supported.
- POST - MID: median +0.023, 6/11 positive, SD 0.046, exact p 0.966, Holm 1.0 -> NOT supported.
A3 edge flicker during rest, dark - light rate (events per hour of flicker-tolerant rest):
- median -0.002, 5/11 positive, SD 0.99, exact p 0.966 -> NOT supported: no day/night difference, so no
  evidence that flicker during rest is behaviour (consistent with a beam-boundary effect; a null at n = 11 is
  weak evidence). Secondary, middle zone only: median +0.017, 8/11 positive, p 0.11.
  Secondary, sleep-measure sensitivity: not counting flickers as movement adds 1-5 min/day of sleep in 10/11
  flies (36 min/day in ch21) and changes bout counts by <= 3% except ch21 (+7%): the 5-min sleep measure is
  robust to flicker in this pilot.
A4 walking, dark - light (Holm over 2):
- median walking speed: median difference -0.231 beams/s, 0/10 positive, SD 0.20, p 0.0077 (normal
  approximation, ties), Holm 0.015 -> SUPPORTED: flies walk more slowly at night than by day.
- median path length per walking bout: median difference 0.0 beams, 3/10 positive, Holm 0.10 -> not supported.
Week 1 (descriptive, code 3fb5405, same measures): A2 PRE - MID +0.092 with 10/10 positive (exact p 0.002),
POST - MID +0.045 (7/10); A3 median +0.013 (6/12); A4 speed -0.105 (1/12 positive, p 0.074), path 0.0.
NOTE: week 1's strong A2 wind-down pattern did NOT replicate in weeks 2-3 under the pre-specified test; A4's
direction did. A descriptive pattern in viewed data is not a finding; this is the case in point.
Sleep per day in these flies (weeks 2-3, standard measure): 726-1229 min.

## 2026-11-14 — A5 ageing (code 2b9ac4f; PRESPEC A5 + P1-F; week-6 phase by hand check)
PRIMARY TEST, within-fly slopes over weeks 2-5 (9 flies with >= 3 qualifying weeks), Wilcoxon two-sided, Holm over 4:
- M1 sleep min/day: median slope +28.0 min/day per week, 7/9 positive, SD 24.8, exact p 0.0195, Holm 0.078
- M2 changes per waking hour: median +9.0, 7/9 positive, p 0.43, Holm 0.43
- M3 light-phase D (A1): median +0.031 per week, 8/9 positive, p 0.0195, Holm 0.078
- M4 median walking speed: median +0.073 beams/s per week, 7/9 positive, p 0.035 (approx.), Holm 0.078
-> NONE supported after Holm (three nominal p < 0.05, all in the same direction: more sleep, more end-resting by
   day, faster walking as flies age from ~1 to ~5 weeks). With 9 flies the smallest exact Holm-adjusted p
   attainable for the first-ranked measure is 4 x 2/512 = 0.016.
P1-F deaths vs the survival sheet: ch17 wk5 (end of file), ch18 wk7, ch20 wk8, ch21 wk6, ch22 wk7, ch23 wk8, ch24 wk8,
ch26 wk6, ch28 wk6, ch31 wk7; undetermined: ch25 (after wk4; wk5 unreadable), ch29 (after wk8), ch32 (from wk1);
never dead: ch19 (alive in wk9), ch27 and ch30 (unreadable, as the sheet implies). One disagreement: ch25, sheet
"died in/after wk5", rule undetermined (its week 5 is unreadable). The ch27 mislabel of P1-E is gone.
DESCRIPTIVE (survivor-biased; flies per week 12, 4, 11, 8, 10, 5, 4, 1, 2): mean sleep 1009 / 922 / 1016 / 1050 /
1083 / 1079 / 1102 min/day in weeks 1-7; light-phase D 0.28 -> 0.49 (wk6); walking speed 1.23 -> 2.03 beams/s (wk6).
Death-aligned (last 7 days vs days 8-21 before death, 9 flies, few valid days near death because of the first-day
and transfer exclusions): sleep rose in 6/8 flies with data, activity per waking hour fell in 5/8.
Week 2 holds only 4 flies: five channels read nothing that week (whole-week outages).

## 2026-11-14 — A2b PRIMARY TEST, weeks 2-3 (code c46431d; PRESPEC A2b = port of hmm_audit PRESPEC 21 + amendment B, LD-adapted)
11 flies (none excluded for < 5 bouts); bouts used per fly 14-147 (L300), 46-262 (L120); only bouts >= 1 h from a light switch.
PRIMARY setting L = 300 s, 10-s bins:
- (a) per-bout null: z median +0.37, Wilcoxon p 0.042, flies with p_fly < 0.05: 2/11 (binomial p 0.10) -> FAILS (binomial arm)
- (b) per-fly day-block: z median +0.29, Wilcoxon p 0.12; (c) group fly-day block p 0.009
-> VERDICT: "no shape asymmetry detected".
SECONDARY setting L = 120 s, 5-s bins (reported; by the rule it does not change the verdict):
- (a) z median +2.15, Wilcoxon p 0.002, 6/11 flies p_fly < 0.05 (binomial p 6e-6); (b) z median +1.72, Wilcoxon p 0.002;
  (c) group block p 0.0005 -> all three arms pass ("TIME-ASYMMETRIC SHAPE" in the secondary setting). Count-intensity
  shape agrees (z median +1.93, p 0.002). Level (first 60 s, start - end): median -0.014, p 0.32.
DIRECTION (descriptive, median per-fly d(tau) = activity after waking - activity before sleep onset, 5-s bins from the
boundary, L120): -0.028, -0.058, -0.026, -0.011, -0.007, then +0.005 to +0.019 from ~25 s on. I.e. in the first ~25 s
after waking flies are LESS active than in the last ~25 s before falling asleep, and more active beyond that. Read
plainly: waking up builds up over ~20-30 s, while activity continues until shortly before sleep onset. This is the
OPPOSITE of the earlier project's description (per0 DD / CS DD, DAM 1-s counts: gradual wind-down into sleep over ~5 min,
abrupt return within ~20 s). Week 1 (descriptive) shows the same near-boundary sign (-0.027, -0.037, -0.024).
Caveats: secondary setting only; LD (light and clock drive not excluded); different detector (beam crossings vs counts),
genotype context (CS LD here) and data from a different experiment.

## 2026-11-14 — A2b port reproduces the earlier per0 DD wind-down result exactly (beam_hmm/repro_a2b_per0dd.py)
Data: the 25 per0 DD 1-s files (per0DD_activity_files_1sec_res, uploaded by the user; not in the repo), cut to the
trace lengths in analysis_20r/per0DD_bouts_20r_replication.mat. Both the ORIGINAL planB_21.py functions and the A2b
port (pipeline.sleep_bouts + the same windows and statistics) give, identically, and identical to planB_21.log:
L300 z median +0.58, Wilcoxon 0.0136, flies p<.05 0.32, day-block z median +0.93 (p 0.0051), d(tau) +0.002 +0.030 +0.045 ...;
L120 z median +0.95, Wilcoxon 0.000912, flies p<.05 0.44, day-block +0.66 (p 0.0275), d(tau) -0.007 +0.019 +0.019 ...
So A2b's opposite-direction profile on the pilot is not a porting error. The earlier per0 DD files are ALSO multibeam
recordings (file names "multibeam_pruned"), so "beam vs count detector" does not explain the difference either.
Remaining candidate explanations (untested): the preprocessing that turns 0.4-s positions into 1-s activity (the pilot
uses spike filter + hub repair + "any change in the second"; the per0 DD files come from the lab pipeline, PreProcess.m
style), LD vs DD, genotype/cohort, and fly age.

## 2026-11-14 — POST HOC: A2b with the LAB activity conversion (code a4d4c8d; labconv.py = port of PreProcess.m +
## multibeam_locomotion.m, binner 5; same 11 flies, valid seconds, ZT and statistics as the A2b test, weeks 2-3)
Not pre-specified; asked to separate preprocessing from biology after A2b's direction disagreed with the earlier project.
- L300: z median +0.61, Wilcoxon 0.014, flies p<.05 2/11 (binomial 0.10, fails), day-block z +0.42 (p 0.032).
- L120: z median +2.56, Wilcoxon 0.00098, 6/11 flies p<.05, day-block z +1.65 (p 0.00098).
- Profile d(tau) (L120, 5-s bins): -0.009, -0.057, -0.025, -0.004, -0.002, then +0.005 to +0.021: the SAME direction as with
  the beam_hmm conversion (less active in the first ~20-25 s after waking than before falling asleep).
=> The lab conversion does not reverse the pilot's direction. Preprocessing does not explain the difference from the
   earlier per0 DD / CS DD result; LD vs DD, genotype, cohort or age remain.
Side result for the lab pipeline (sleep per day, same seconds): the two conversions differ in 0.3-5.7% of seconds; sleep
agrees within ~1% in most fly-weeks, but the lab conversion gives much LESS sleep where a blocked-beam (hub) artifact is
present, because it counts the artifact as movement: ch20 wk2 966 vs 1227 min/day, ch29 wk3 595 vs 763, ch17 wk2 1001 vs
1041 (and ch21 wk2 679 vs 621 the other way, where the beyond-beam-1 zone counts as movement only in beam_hmm's version).
Note: labconv.py is checked against the MATLAB code by reading only; a MATLAB golden run on one raw file would confirm it.

## 2026-11-14 — Lab sleep scoring vs beam_hmm scoring: identical (labconv.sleeptrace/ironout = ports of sleeptrace.m,
## ironout.m; process_LocomToSleep.m read)
sleeptrace (5-min rule) + ironout (1-s interior wake runs flanked by sleep -> sleep, one pass) vs pipeline.sleep_bouts
(sqw_ref rule: a single active second is bridged iff the unbridged inactive runs on both sides are each >= 300 s):
weeks 2-3, every live channel 17-32, both activity conversions (beam_hmm and lab) -> 50 fly-week-conversions,
27,978,555 seconds, 0 differ. So all differences between the lab pipeline and beam_hmm lie OUTSIDE sleep scoring:
1 activity conversion (PreProcess + multibeam_locomotion vs spike filter + hub repair + 0 kept): measured above;
  sleep differs by up to 261 min/day in fly-weeks with a blocked-beam (hub) artifact;
2 death rule: process_LocomToSleep cuts 24 h before the FIRST >= 12-h run of zero activity in each file; beam_hmm (P1-F)
  works across files and separates dead from unreadable. Applied to this run's week 1, the lab rule would treat the
  ~14 h with tubes out of the monitor (from ~153 h) as a death in every fly and cut from ~129 h; and a blocked-beam
  channel never reads zero, so it cannot be called dead;
3 phase: startZT = 0 (file start) in the lab call vs the flies' own lights-off rise here (files 1, 2, 5, 7 start
  ~30-40 min away from lights-on);
4 entrainment: Teliminate 24 h = beam_hmm's first-24 h exclusion (same).

## 2026-11-14 — POST HOC: A2b without the light-switch exclusion (code 56f4013, --no-switch-exclusion; weeks 2-3, 11 flies)
All wake bouts used (as in the earlier project's CS LD block). Bouts per fly 25-190 (L300), 68-320 (L120).
- L300: z median +0.20, Wilcoxon 0.52, flies p<.05 2/11 -> no shape asymmetry (as with the exclusion).
- L120: z median +3.13, Wilcoxon 0.00098, 7/11 flies p<.05, day-block +2.57 (p 0.00098), group block p 0.0005 -> all arms pass.
- Profile d(tau) (L120): -0.028, -0.056, -0.036, -0.019, -0.009, then +0.009 to +0.021: the SAME direction as with the exclusion.
=> The light-switch exclusion does not explain the pilot's direction. Every method difference with the earlier project is
   now ruled out; the remaining explanations are biological or experimental (cohort, age, genotype, LD vs DD).

## 2026-11-14 — Food-blocking check of hub beams; RawToActivity1s golden test
Hub beams per fly-week (spike-filtered, weeks 1-9, channels 17-32): 35 fly-weeks where the dead-hub rule fired;
13 at an end beam (beam 16 in ch17, 20, 22, 27, 30, 31, 32; beam 1 in ch26), 22 mid-tube (beams 3-11). End-beam hubs
often read 60-100% of samples (fly nearly invisible). Within the week, the hub's share of readings fell from day 1 to
day 6 in 8/13 end-beam cases (mean 0.61 -> 0.56; e.g. ch17 wk2 0.45 -> 0.03, ch30 wk3 0.83 -> 0.62, ch32 wk6 0.99 ->
0.60): consistent with food (or condensation) blocking an end beam and clearing, NOT decisive (descriptive, post hoc;
food end per tube unknown). Mid-tube hubs cannot be food at an end.
Runs of 0 (all weeks): 109,618, every one entered from beam 1, none from beam 16: "no beam broken" occurs only beyond
beam 1. Open: the monitor's rule when two beams are broken, and the food end per tube.
RawToActivity1s.m (beam_hmm/matlab): golden test vs labconv.py passes for all 4 variants x 16 channels on the one-day
excerpt in Octave 8.4; full-week MATLAB run pending on the user's machine.

## 2026-11-14 — LoadingCheck (first 2 h after loading) validated against the full-week hub QC (weeks 1-9, ch17-32)
LoadingCheck.m and its Python twin (labconv.loading_check) give identical flags on 64 channel x 2-h excerpts (Octave 8.4).
BLOCKED BEAM flag in the first 2 h vs "hub declared dead" over the whole week, 144 fly-weeks: both 18; 2-h only 1
(wk1 ch24, beam 4); week only 20; neither 105. Of the 20 missed, 9 were flagged LOW ACTIVITY in the first 2 h (so a
re-seat check would still have caught the tube); 11 had no flag -- hubs that develop or strengthen later in the week
(e.g. wk1 ch17: hub-16 readings 0.00 on day 1, 0.04-0.07 later) or give too few jumps in 2 h.
Reading: the 2-h screen catches ~half of blocked beams directly and ~70% with any flag; it complements, not replaces,
the full-file QC (RawToActivity1s / BeamFingerprint). Thresholds were NOT tuned on this run (they are the P1-C hub
criteria plus a minimum of 10 jumps); tuning them would need a new, dated entry and new data.
In weeks 8-9 nearly every channel is LOW ACTIVITY: the flies are dead, as the survival record says.
    CORRECTION (same day): the "64 channel x 2-h excerpts" used for the Octave-vs-Python equality check above held
    only ~1 h (the line count ignored the alternating -1 lines), and LoadingCheck.m's own line budget was too small
    in Octave (it counts CRLF line ends twice). Fixed (budget doubled; the time window still trims to 2 h); rerun on
    true 2-h excerpts of weeks 1-4: identical flags, 0 short-file warnings. The validation against the full-week QC
    (18 / 1 / 20 / 105) used the full data and is unaffected.

## 2026-11-14 — Live flies hidden behind a blocked beam 16; LoadingCheck FROZEN flag
Channel-weeks with one value in >= 99% of readings (weeks 1-9, ch17-32): 44 -- value 16 in 25, value 0 in 17 (dead flies,
e.g. all of week 9), 2 and 6 once each. In 14 of them the reading was a constant 16 for the WHOLE week while the fly
was alive (active again in a later week): wk1 ch30; wk2 ch19, 27, 29, 30, 31, 32; ch27 in weeks 3-8; ch32 wk5. A blocked
beam 16 (plug or food, the user's hypothesis) can hide a living fly completely: ~1 in 10 fly-weeks of the pilot.
These are the "whole-week outages" of A5; P1-F already treats them as outages, not deaths.
LoadingCheck/labconv.loading_check gain FROZEN AT END BEAM k / FROZEN AT BEAM k / FROZEN AT 0. In the first 2 h it
flags all 14 of 14 hidden-alive fly-weeks as FROZEN AT END BEAM 16. Octave and Python flags identical on 2-h excerpts.

## 2026-11-14 — T1-B synthetic gate (code dc3a50c + T1-B1; run on the user's machine, pushed as f5eb874): PASS
12-fly cohorts x 14 days, per-fly heterogeneity, K fixed at the true 3, sample-level HMM.
B-P (planted morning/evening sleep-exit rhythm, peak-to-trough 3): C2 vs C0 rejected in 20/20 (need >= 16; chi-square
p 1e-236 to 3e-178, LR median 961); dark-phase peak of the fitted exit profile within 2 h of ZT 12.5 in 20/20 (ZT 12.3-12.7);
fitted/true within-dark range median 1.03 (0.92-1.11; need 0.67-1.5). B-N (no rhythm): 0/20 rejected (allowed <= 2;
chi-square p mean 0.37, 2/20 below 0.10). Sleep state identified as QUIESCENT/S in all 480 fly fits (0 fallbacks).
Bootstrap vs chi-square (5 + 5 cohorts, 99 replicates): B-P p_boot 0.010 (the minimum attainable) in 5/5; B-N p_boot
0.25-0.69 vs chi-square 0.14-0.71 -- chi-square p was the smaller in 4/5 null cohorts (mildly anti-conservative), so the
pre-specified real-data decision by the bootstrap stands.
Limits, stated: the gate fixes K at the truth and the simulator's dynamics are close to the fitted model family; it
does not test state discovery (D1-D failed; B1 pending). Power at the planted 3-fold rhythm is saturated, so the gate
says nothing about smaller rhythms. Per T1-B, a real-data run needs a state model that has passed its own gate.

## 2026-11-14 — T1-C synthetic gate (code b4dd91d; run here): FAIL (false-detection arm)
B-P (planted 3-fold sleep-exit rhythm): rejected in 20/20 (rotation p 0.001, the minimum, in all); dark-phase peak of the
fitted wake-up shape within 2 h of ZT 12.5 in 20/20 (ZT 12.2-12.8) -> PASS.
B-N (no rhythm beyond the light/dark step): rejected in 3/20 (p 0.020, 0.036, 0.041; allowed <= 2) -> FAIL.
Null p values: mean 0.46, 5/20 below 0.10. 3/20 is not strong evidence of miscalibration on its own (P(X >= 3) = 0.075
at a true size of 0.05), but the pre-specified rule fails and stands.
Likely mechanism (stated, not tested): in the observed scoring, the wake-up hazard depends on time asleep (bouts last
>= 5 min), so after each light switch the mix of short and long sleepers re-equilibrates; that transient is aligned to
ZT in every fly, i.e. SHARED within-phase structure produced by the step itself, which the rotation null destroys.
The hidden-state T1-B model has no such transient (its rates are constant within a phase), which is why T1-B passed.

## 2026-11-14 — T1-C1 implementation note (before any T1-C1 gate result was read)
The first T1-C1 gate start was stopped after 4/80 cohorts (no result read) and restarted from scratch with the same seeds:
(1) the Poisson fits now maximise the PROFILE likelihood over the shared harmonics, with the per-fly x phase intercepts
in closed form -- exactly the same maximum; on a simulated cohort it agreed with the full fit to <= 8.4e-4 in the LR
(the full L-BFGS fit stops slightly short; the profile fit is the tighter one); (2) worker processes run with
OMP/OPENBLAS/MKL threads = 1 (4 workers x 4 BLAS threads on 4 cores had slowed the run ~50x).

## 2026-11-14 — T1-C1 synthetic gate (code 5463e67; run here): PASS
Fresh cohorts (seeds 9501 / 11201), 12 flies x 12 days. E = 2 h (decision): B-P rejected 20/20, dark-phase wake-up peak
within 2 h of ZT 12.5 in 20/20; B-N rejected 2/60 (0.033; allowed <= 5). Null p quartiles 0.28 / 0.61 / 0.85.
Sensitivity (same data): E = 1 h B-N 2/60, E = 3 h 3/60; B-P 20/20 at both. The exclusion window barely matters on
these cohorts, so T1-C's 3/20 was probably partly chance as well as the switch transient; T1-C1 stands as specified.

## 2026-11-14 — T1-C1 REAL-DATA TEST, weeks 2-3 (code 3472418; PRESPEC T1-C + T1-C1): SUPPORTED
11 flies (ch17-24, 26, 28, 29; none dies inside the window). E = 2 h (decision): LR 1012.5, rotation p 0.001 (the minimum
with 999 rotations; null median LR 81.4) -> time-of-day dependence of sleep onsets and wake-ups beyond the light/dark
step, shared by the cohort, is SUPPORTED (LD: light + clock). Sensitivity: E = 1 h LR 1002 (p 0.001), E = 3 h LR 977
(p 0.001). Secondary C1 (24-h harmonic only) vs C0: LR 177.5, p 0.019.
Fitted within-phase shapes (E = 2 h; descriptive): wake-ups peak late in the light phase (ZT 10.8) and at the start of the
dark phase (ZT 12.1 -- inside the excluded 2 h, i.e. the harmonic fit's extrapolation toward lights-off); onsets peak at
midday (ZT 5.1, the siesta) and in the early night (ZT 16.9). Dark-phase log-rate range 0.80 (wake-ups), 0.98 (onsets).
Week 1 (descriptive): same picture (LR 487, p 0.001; wake-up peak ZT 11.9 / 12.1; onset peaks ZT 5.2 / 17.2).

## 2026-11-15 — PLANNING (not a result): sample sizes for the repeat experiment (code 6633f48 design_power.py)
Uses only per-fly summaries already in this ledger (pilot median / fly-to-fly SD) and simulates the pre-specified
one-sample Wilcoxon at the Holm worst case (alpha/m), 4000 replicates per n. Analysable flies for 80% power at the
pilot effect / at half of it: A1 light D 8 / 15; A4 speed 12 / 40; A5 M1 15 / 40; A5 M3 50 / >100; A5 M4 100 / >100;
A5 M2 >100 / >100; A2 illustrative 0.03 effect 60 / >100. Output results/design/power.json; reading in DESIGN_NOTE.md.
No new data were looked at.

## 2026-11-15 — LoadingCheck 'StartHour' / 'Previous' added; first-screen window 1 h vs 2 h (code 3fc1e50, loadingcheck_window.py)
New options: 'StartHour' (window start after the file's first reading, to re-check after a re-seat in a running file) and
'Previous' (an earlier R; per channel RESOLVED / STILL FLAGGED / NEW FLAG / OK). Octave LoadingCheck.m vs Python twin on
weeks 1-2, ch17-32 (first 2 h, then a 1-h window at 30 h with 'Previous'): identical flags and labels, 32/32.
Window comparison, weeks 1-9, ch17-32 (144 fly-weeks; thresholds unchanged):
- 2 h reproduces the 2026-11-14 validation exactly: BLOCKED vs full-week dead hub 18 both / 1 window only / 20 week only
  (9 of them with another flag) / 105 neither; hidden live flies 14/14 FROZEN AT END BEAM 16.
- 1 h: identical on all of these.
- Flags on READABLE fly-weeks with no full-week hub and no hidden fly (59): 7 at 2 h (LOW ACTIVITY 6, FROZEN 4, BLOCKED 1),
  8 at 1 h (+1 NARROW RANGE). Candidate false alarms, not proven ones.
Side observation (the 30-h Octave check, week 1 = ZT ~5, siesta): several active flies came out LOW ACTIVITY and one
FROZEN AT END BEAM 16 -- a window during sleep flags sleeping flies. Re-checks should follow the re-seat promptly.
Recommendation unchanged: 2 h (default); 1 h acceptable for a first look on this evidence (one pilot run).

## 2026-11-16 — LoadingCheck re-check options verified in MATLAB R2026a (user's PC at 40de248; PNGs pushed as f3ede2c)
Week 2, ch17-32: R1 = first 2 h; R2 = 2.5-4.5 h with 'Previous', R1; R3 = 1-h window. R2 flags and Change labels match the
Python twin on all 16 channels (ch17 RESOLVED; ch24 NEW FLAG; 9 STILL FLAGGED; 5 OK). Printed summary "1 resolved, 9 still
flagged, 1 new, 5 OK" -- correct; the expected line sent from here ("10 still flagged, 4 OK") was my counting slip. R3 printed
the ShortWindow warning; no errors; MATLAB -batch exited normally. Figure: white background, RESOLVED green, STILL/NEW red,
x label "min after 2.50 h".
Not a re-seat experiment: nothing was re-seated at 2.3 h in the pilot, so this tests the code only. ch24's NEW FLAG (FROZEN AT
END BEAM 16; LOW ACTIVITY) shows the trace sitting at 16 with one excursion: a fly resting at the end beam and a fly hidden
behind a blocked end beam look the same in a window of quiet behaviour.
Later commit b02159e (made from the user's Box analysis folder by another session: 16-channel figure pages, file dialog,
full-path output, 'Previous' re-reads its file, review fixes) leaves the flag logic unchanged: re-checked here in Octave
against the Python twin on weeks 1-2 (2 h, and 1 h at 30 h with 'Previous'), 32/32 identical.

## 2026-11-16 — FACT (user): position 0 = at the food in a correctly oriented tube; DESCRIPTIVE look at existing A1 output
User: with the tube oriented correctly, beam 1 is the beam nearest the food, so 0 means the fly is very likely eating.
Open: the info sheet says the food end differs between flies (PROTOCOL.md) -- were some pilot tubes reversed?
Descriptive only (numbers already in results/a1/a1_test_weeks2-3.json, the A1 BEYOND1 secondary; nothing recomputed),
12 primary flies, weeks 2-3, share of scored SLEEP / PAUSE samples at position 0:
- light: sleep median 0.000 (max 0.029, ch21); pauses median 0.010 (max 0.040)
- dark:  sleep median 0.002 (max 0.075, ch21); pauses median 0.010 (max 0.361, ch21)
Reading: little scored sleep happens at 0, so "a still fly at the food scored as asleep" moves the pilot's sleep totals
little (at most ~7.5% of one fly's dark sleep). Short stillness at 0 is ~1% of pauses (36% of ch21's dark pauses).
