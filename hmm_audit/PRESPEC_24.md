# PRESPEC 24 — 23E10-GAL4 silencing: active Kir2.1 against non-conducting Kir2.1 (same background)

Date 2026-11-06. Written and committed BEFORE any behavioural quantity of these flies has been computed. Only file
metadata (counts, lengths, names) has been inspected. This is an exception to the freeze (new experiment).

## Data (Sheyum, uploads of 2026-11-06)
- 23E10_x_Kir.zip (SHA-256 33449da5...): 59 flies, 23E10-GAL4 > UAS-Kir2.1 (active), monitor "4brc2", 750,024 s each.
- 23E10_x_NC.zip (SHA-256 6ee4d9fe...): 64 flies, 23E10-GAL4 > UAS-Kir2.1 non-conducting (control), monitor "4bWC",
  1,082,745 s each.
- Same genetic background apart from the active/mutated Kir. Files start at ZT0. Assumed (as for 20y, to be
  confirmed by Sheyum): males, LD 12:12. Different monitors: monitor and genotype are confounded, as stated with every
  result. Raw data are not uploaded to the repository.

## Records (fixed now)
- Death (Sheyum's rule): the first run of >= 12 h with no activity marks death; the record is cut 24 h before its
  start. Whatever remains is used.
- Days: the first 8 whole days (the length of the Kir records), applied to both groups after the death cut. A fly with
  fewer than 2 whole days left is excluded. Exclusions are counted per group.
- Symbols and readouts: identical to 20y (sqw20y.one): sqw_symbols theta 0.5, light/dark references, 12-h phases,
  observable_quantities (c_obs, count J, count EPR), sleep min/h, the six transition probabilities, E_W - E_Q and
  E_S - E_Q.

## PRIMARY (confirmatory; one test)
From 20y, dFB split-GAL4 silencing raised dark-phase J by +0.104 cycles/h (MW p 0.055, not significant after Holm). 23E10
targets dFB neurons, and reportedly also sleep-promoting VNC neurons. Hypothesis: silencing raises dark-phase cycle
flux. Test: one-sided Mann-Whitney, 23E10>Kir > 23E10>NC, on count J in the dark phase (flies with a dark
reference). Alpha 0.05, no correction. Effect: median difference with fly-bootstrap 95% CI, plus Cliff's delta.

## SECONDARY (two-sided)
- Family A (Holm over 6): c_obs, J, EPR x light, dark. J dark appears here too, two-sided.
- Sleep min/h, light and dark (Holm over 2).
- Six transition probabilities x 2 phases (Holm over 12); E_W - E_Q and E_S - E_Q x 2 phases, cycle flies (Holm over 4).
- Sensitivity: everything without day 0 (as the 20y 'nod0').

## Across backgrounds (descriptive only; no pooled test)
For each readout, the sign and size of the median difference against the group's own control, side by side for 23E10
(this experiment), dFB split and VNC split (20y, against empty>Kir). Reported as concordant or discordant in sign. Values
are never compared across backgrounds.

## Reading (fixed now)
- Primary p < 0.05: "silencing 23E10 neurons raises the night-time cycle flux". This confirms the 20y dFB hint in an
  independent, background-matched experiment.
- Otherwise: "not confirmed". The 20y hint is then reported as not replicated.
- Any sleep change is interpreted with C14 in mind: less sleep alone LOWERS J (fewer transitions). A J increase
  together with less sleep therefore points to a change of the drive, not of sleep amount.

## Output
`sqw24.py` -> `sqw_out/s24_23E10.json`, `figures/fig24_23E10.png`. Ledger entry with commit hashes.
