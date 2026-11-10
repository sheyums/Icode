# PRESPEC 30 — does sleep deprivation increase LONG-bout sleep? (Chowdhury et al. 2023 criterion)

Date 2026-11-10. This is an exception to the freeze, approved by Sheyum. Written and committed BEFORE the deprivation
recordings were available in this session (they reached cClaudeII for 20x; nothing in this session has computed any
bout-level quantity from them). Already known from 20x: no change of the cycle or of all-bout sleep in ZT0-6 beyond
the day-to-day drift; ZT0-6 sleep rises through the whole recording.

## Motivation
Chowdhury, Abhilash, Ortega, Liu & Shafer (2023) eLife 12:e91355: homeostatic rebound is carried by long sleep bouts
(inactivity criterion >= 25 min gives the largest payback, about 50% by recovery day 3). Bout DURATION rises after
deprivation, bout number does not. The rebound spans several days and the whole day. It was detected only against
yoked controls (same mechanical stimulation), not against undisturbed flies.

## Data (as 20x)
CS_9hrdeprExp: 26 CS males, 9.513 days each, from ZT0. Deprivation day index 3, ZT15-24 (excluded from every
readout). Baseline B = days 0-2; recovery days 4, 5, 6. Placebo (the drift control; there are no undeprived or yoked
flies): day 2 against days 0-1.

## Bouts
Sleep bouts as sqw_ref.sleep_bouts (bridged 5-min criterion, 1-s resolution). LONG bout = a sleep bout lasting
>= 25 min (primary) or >= 30 min (fixed second criterion, reported). Long-bout sleep in a window = the seconds of
long bouts falling in that window (a bout crossing midnight is split by second). A bout's duration, start and exit
belong to the day of its start.

## Per-fly readouts and contrasts
For readout X: Delta = X(recovery day) - mean X(days 0-2); Delta_pl = X(day 2) - mean X(days 0-1).
PRIMARY: X = long-bout (>= 25 min) sleep in minutes over the whole day (ZT0-24), recovery day 4.
  A rebound is claimed only if BOTH Wilcoxon(Delta) > 0 (p < 0.05) AND the paired Wilcoxon of Delta - Delta_pl > 0
  (p < 0.05). Two-sided tests; the sign must be positive.
SECONDARY (the Delta - Delta_pl test; Holm over 8):
  S1, S2 long-bout sleep on recovery days 5 and 6;
  S3 the >= 30 min criterion, day 4;
  S4 the long-bout share of sleep (long-bout sleep / all sleep), day 4;
  S5 the mean duration of sleep bouts (all bouts >= 5 min) starting on day 4;
  S6 long-bout sleep in ZT0-6 of day 4 (the 20x window);
  S7 the exit route of long bouts on day 4: the fraction whose first wake minute is active wake (W), on 20x's
     baseline-referenced symbols;
  S8 the entry route of long bouts: the fraction whose last wake minute before onset is quiet wake (Q).
Descriptive: per-day time course of long-bout and all-bout sleep, days 0-8; the number of long bouts per day.

## Reading (fixed now)
- Primary pass: "long-bout sleep rises after deprivation beyond the day-to-day drift". This is consistent with a
  homeostatic response confined to long bouts, but mechanical side effects are not controlled (no yoked flies).
- Primary fail: "no long-bout rebound detectable beyond the drift". Following Chowdhury et al., this is WEAK evidence:
  their rebound was invisible without yoked controls. The controlled experiment should pre-register this readout,
  with yoked controls.
- S7/S8 positive: deprivation shifts how long bouts are entered or left (the cycle's forward step, C7). Hypothesis
  for the controlled experiment, not a claim.
- No other bout threshold, day, window or readout is added after looking.

## Output
`sqw30.py` -> `sqw_out/s30_longbouts.json`, `figures/fig30_longbouts.png`. Ledger entry with commit hashes.
