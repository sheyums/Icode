# PRESPEC 20x — merged paper, test 4: does sleep pressure change the observed S/Q/W cycle? (before running)
Date 2026-11-01. Approved by Sheyum for the merged paper; brief by cClaude. Owner cClaudeII, branch planB-feasibility.
Results go to RESULTS_LEDGER_planB.md. The F4 check (20t) has already looked at the group hourly sleep profile and
per-fly sleep lost; nothing else from this dataset has been computed.

## Question
Does sleep pressure change the observed cycle (its direction, flux and dissipation), and which transitions does it
change? Deprivation is the homeostat's counterpart to the clock result of 20w-B. Per-fly regulator fits failed
feasibility (20t F1), so the design uses within-fly contrasts summarised at group level.

## Data
CS_9hrdeprExp.zip: 26 flies, 9.513 d each, from ZT0. Deprivation on day index 3, ZT15-24 (located in 20t F4 from the
group profile; it matches Sheyum's "last 9 h of the night"). Days 0-2 are baseline; days 4-8 recovery. The deprivation
night (day 3 ZT15-24) is excluded from every readout, because mechanical arousal makes its symbols meaningless.

## Symbols (one stated change from sqw_symbols)
sqw_ref sleep scoring (bridged 5-min criterion, S if >= 30 s of the minute in a sleep bout), theta 0.5, light/dark
references. The change: each fly's r_f (light and dark) is computed from wake-bout events of its BASELINE days 0-2
only (sqw_ref.wake_bout_events on the first 3 days; wake_reference per light/dark block) and applied to every day.
Otherwise the recovery activity level would shift the Q/W boundary. Implemented as a thin wrapper reproducing
sqw_symbols with fixed references. A wake minute in a block without a reference is missing (-1), as in sqw_symbols.
Reported: the shortest trace, and the number of flies without a reference in each light/dark block.

## Window (fixed before any data)
- PRIMARY: ZT0-6, minutes [0, 360) of the day. Sheyum chose it in the design discussion BEFORE the data were uploaded.
  20t F4 has since seen a group rebound at ZT1-6. Both facts are recorded: the window was not chosen after looking.
- DESCRIPTIVE only: day 4 ZT6-12 and ZT12-24 (each against the same window of days 0-2), and ZT0-6 of days 5-8
  (time course of return to baseline). No other windows; no widening (rule 4).

## Per-fly contrasts (all in ZT0-6)
B = days 0-2 pooled; R = day 4; placebo (time-trend control; there are no undeprived flies) P = day 2 against days
0-1 pooled. P has the same window length and the same flies, but one fewer baseline day.
Delta = R - B; Delta_pl = day2 - days01.
A readout on pooled minutes uses the minutes of all its days together, with transitions attributed to the window of
the arrival minute (sqw_ref convention), and hours = minutes in the window / 60.

## Stage 0 gate: did the deprivation work?
Per fly, Delta of minutes of sleep in ZT0-6 (day 4 minus the mean over days 0-2), sleep = seconds inside bridged
sleep bouts / 60. PASS if the group two-sided Wilcoxon gives median Delta > 0 at p < 0.05. FAIL => report
"no measurable rebound" and run nothing more.

## Primary readouts (family of 3; Holm across them)
c_obs (sqw_ref.c_obs_minutes), count J (sqw_ref.j_cycle_count, cycles/h), count EPR (sqw_ref.epr_count_bound, rule of
three). Count EPR is finite by construction, so no Inf can occur. Medians would rank Inf highest if one did.
"No cycle" fly-conditions are classified into the three kinds of the 20s-A correction and reported:
(i) one-way edge (rule of three applied); (ii) exactly symmetric counts (net J = 0); (iii) an edge absent both ways or
a state absent in the window.
DECISION per readout, SUPPORTED only if BOTH hold in the same direction:
(1) two-sided Wilcoxon of Delta against 0 at Holm p < 0.05;
(2) paired Wilcoxon of Delta against Delta_pl at p < 0.05, with median(Delta - Delta_pl) of the same sign as median Delta.
The Delta_pl test against 0 is reported as the time-trend check; if significant it is FLAGGED.

## Secondary: which transitions
Group-pooled observed chain on the ZT0-6 minutes: transition counts between consecutive non-missing minutes, pooled
over flies, B (days 0-2) vs R (day 4), row-normalised. This is mc_fit's count MLE: identical to mc_fit on complete
minutes; missing minutes are left out rather than marginalised, an implementation simplification stated here.
For each of the 6 off-diagonal transition probabilities: ratio R/B with a fly-bootstrap 95% CI (2000 replicates,
resampling flies, seed 2026), bootstrap p = 2 min(frac(ratio* <= 1), frac(ratio* >= 1)), Holm across the six.
Wording: "sleep pressure changes transition X -> Y". Never "the homeostat acts on X -> Y".

## Secondary: dose association
Spearman rho between sleep lost (F4's per-fly number: mean ZT15-24 sleep over days 0-2 minus day 3 ZT15-24 sleep) and
each primary Delta. Computed in all flies, and in efficacious flies only (deprivation-night sleep <= 108 min, i.e.
awake >= 80% of the 540-min window). Dose is partly a fly property (resistance to arousal; F4 caveat), so the wording
is "associated with", never "dose-response".

## Narrowest block, matched null
Narrowest block: one fly's ZT0-6 on day 4 (360 min). The placebo contrast is the sample-size-matched null: same flies,
same window length, same estimator.

## Wording, fixed now
SUPPORTED: "After 9 h of night-time sleep deprivation, [readout] changed by [median Delta, CI] in the first 6 h of
recovery, beyond the day-to-day variation of the same flies."
Not supported: "No change of [readout] beyond day-to-day variation was detected", with the CI quoted (fly-bootstrap 95%
CI of median Delta, 2000 replicates).
One run, one verdict. Time box 2 h.
