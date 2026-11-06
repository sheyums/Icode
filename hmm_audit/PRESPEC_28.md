# PRESPEC 28 — R60D05 (EPG) silencing: Kir2.1 against non-conducting Kir2.1

Date 2026-11-06. Written and committed BEFORE any behavioural quantity of these flies has been computed. Only file
metadata (counts, names, lengths) has been inspected. This is an exception to the freeze (new experiment).

## Data
epg_silenced.zip (SHA-256 6fddf6c2...): 64 files "R60D05_x_BLKir*", 886,449 s each. epg_control.zip (SHA-256
11893ea8...): 64 files "R60D05_x_BLnc*", 886,483 s each. Per Sheyum: males, LD 12:12, start at ZT0. The control is
the same GAL4 driving non-conducting Kir2.1 (from the file names; to be confirmed). Near-identical lengths suggest one
run. Raw data are not uploaded.

## Records and readouts: identical to PRESPEC 24
- Death rule: 12 h without activity; cut 24 h before.
- First 8 whole days; at least 2, else the fly is excluded.
- 20y readouts: c_obs, count J, count EPR, sleep, six transitions, landscape depths.
- Values are rounded to 10 decimals before ranking.

## PRIMARY
There is no directional hypothesis (Sheyum expects little sleep change), so the primary family is 20y's: c_obs, J, EPR x
light, dark. Two-sided Mann-Whitney, Holm over 6.

## SECONDARY
Sleep (Holm over 2); transitions (Holm over 12); landscape depths (Holm over 4); without-day-0 sensitivity.

## Across backgrounds (descriptive)
The sign of each median difference against its own control, side by side with 23E10 (PRESPEC 24) and the dFB / VNC
split lines (20y).

## Reading
A cycle change is claimed only at Holm < 0.05 in the primary family. A sleep change alone is "sleep only", as in C12.

## Output
`sqw28.py` -> `sqw_out/s28_epg.json`, `figures/fig28_epg_J.png`. Ledger entry with commit hashes.
