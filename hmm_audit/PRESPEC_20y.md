# PRESPEC 20y — circuit silencing: which parts of the observed S/Q/W cycle depend on dFB and VNC neurons?
Date 2026-11-01. Written and committed BEFORE any behavioural quantity of these flies has been computed. Only file
metadata (names, lengths, format, duplicates) has been inspected. Results go to RESULTS_LEDGER.md (hmm_audit).

## Data (Sheyum, uploads of 2026-11-01)
Constitutive UAS-Kir2.1 driven by: an enhancerless GAL4 (empty>Kir, control), a dFB split-GAL4 (dfb>Kir) and a VNC
split-GAL4 (vnc>Kir). All three in ONE common genetic background (not CS), males, LD 12:12, one run, one incubator.
Files: 1-s counts, one per line, starting at ZT0 (Sheyum), names "...cutdata_4.8.2025.txt".
  empty_kir.zip (SHA-256 01d30b98...): 19 flies; monitor rc2 ch 1-3 (715,157 s) and monitor rc ch 5-30 (715,166 s).
    The two files named "...rc21...20251.txt" and "...rc23...20251.txt" are monitor rc ch 21 and 23 (different flies
    from rc2 ch 1 and 3: 79,048 and 106,856 differing seconds), renamed to avoid a name clash.
  dfb_kir.zip   (SHA-256 d455d42d...): 25 flies; monitor rc2 ch 6-48 (715,157 s).
  vnc_kir.zip   (SHA-256 927f6c7d...): 34 flies; monitor rc2 ch 49-62 (10 flies, 715,157 s) and monitor wc (24 flies,
    715,185 s).
Not uploaded to the repository (raw data rule).

## What can and cannot be compared (fixed now)
- dfb>Kir and vnc>Kir are compared ONLY with empty>Kir. Nothing is compared with the CS / mutant cohorts (different
  background).
- Kir2.1 is constitutive: effects may arise in development. Wording: "constitutive silencing ... changes", never
  "these neurons acutely control".
- The control carries a full enhancerless GAL4; the drivers are split-GAL4s; no empty split-GAL4 control exists. The
  insertions of the split hemidrivers are not controlled. Stated in every result.
- The two papers describing the lines could not be read from this session (network); every test is TWO-SIDED.

## Symbols and records (identical to 20s-A; verified identical on a CS file before writing this)
sqw_ref.sqw_symbols(x, theta=0.5, light_dark=True, full_days=True, unknown_wake=-1): 8 whole days per fly.
Phases: 12 h (ZT 0-12 light, 12-24 dark), by arrival minute. Chain: sqw_ref.mc_fit (2 phases).
Readouts per fly x phase: sqw_ref.observable_quantities (c_obs, count J in cycles/h, count EPR rule of three, chain
J / A / EPR, E21, E31, F_SW, F_WS, Xn, Yn); the six off-diagonal transition probabilities of the fitted chain; sleep
(minutes of S per hour of the phase).
Exclusions, applied mechanically in the run: (a) a fly whose final 24 h contain no active second is called dead and
excluded entirely (listed); (b) a fly x phase without a Q/W reference is excluded from that phase (as 20s-A);
(c) fly x phases with no net cycle current (|J_chain| < 1e-12) are excluded from J-, EPR- and L-based summaries and
counted by kind (symmetric counts / transition never observed / state absent), as 20s-A.

## Primary family (Holm over 12 tests)
Readouts c_obs, count J, count EPR x phases light, dark x contrasts dfb vs empty, vnc vs empty.
Test: two-sided Mann-Whitney U on per-fly values (EPR: Inf ranks highest, as 20s-A medians); exact or normal
approximation with tie correction as scipy chooses. Effect size: probability of superiority U / (n1 n2) and the
difference of group medians with a fly-bootstrap 95% CI (2000 resamples within each group; Inf ranked highest).
Sample-size matching (rule 5): the Mann-Whitney null is a label permutation of the same flies, matched by construction.
Reading: a readout x phase x contrast CHANGES if Holm p < 0.05. Otherwise "no change detected", with the CI.

## Positive control (reported first, NOT a gate)
Sleep per phase, dfb vs empty and vnc vs empty, Mann-Whitney, Holm over 4. A cycle change without a sleep change is
itself a result, so no stopping rule hangs on it.

## Secondary (Holm within each family; descriptive wording "consistent with")
S1 the six transition probabilities x 2 phases x 2 contrasts (Holm over 24).
S2 landscape E21, E31 x 2 phases x 2 contrasts (Holm over 8).
S3 Onsager L per group (EPR = J^2/L through the origin over that group's fly x phases with finite count EPR and a
   cycle, light and dark pooled, as Fig A2), fly-bootstrap 95% CI; difference to empty with a fly-bootstrap p (2000).

## Monitor check (pre-specified, descriptive)
Same genotype, different monitor: vnc>Kir on rc2 (10) vs on wc (24), and empty>Kir on rc2 (3) vs on rc (16), for
sleep, c_obs, J and EPR per phase (Mann-Whitney, no correction). If a monitor difference within a genotype is as large
as a reported genotype difference (by median), that genotype result carries "possibly monitor-related". dfb>Kir is
entirely on rc2; its comparison with empty>Kir rc2 has n = 3 and is NOT used as a test.

## Sensitivity (pre-specified)
Day 0 excluded (acclimation): the primary family recomputed on days 1-7. Reported beside the primary result; the
verdict is the primary one.

## Narrowest block
One fly's light (or dark) minutes pooled over 8 days. Group statistic over flies.

## Bounds
One run, one verdict, no re-binning, no further windows. Time box 2 h.
