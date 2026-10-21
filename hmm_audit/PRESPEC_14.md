# PRESPEC 14 — partial pooling: SHARED routing, FLY-SPECIFIC texture (plan + build + unit test only)
Date 2026-10-21. Motivated by PRESPEC 13: full pooling recovers routing precisely (V1, V5) but fly
heterogeneity creates a false transitional state (V2, V6). By user decision, this step only plans,
builds and unit-tests the model; no validation run and no data run until the design review (pause).

## Model (hhmm_fly.HHMMFly), per genotype x 4 h ZT bin
Coarse c = 1..Nc, fine f = 1..2. For fly i on day d:
  emissions            E_i[(c,f), symbol]              FLY-SPECIFIC
  coarse exit prob     eps_i[c,f] * kappa[c,d]          eps FLY-SPECIFIC, day factor kappa SHARED
  fine dynamics        A_c                              SHARED
  coarse jump matrix   G   (routing: Xn_c, Yn_c)       SHARED
  entry distribution   Bin_c                            SHARED
Fitted by generalized EM: exact closed-form updates for E_i, A_c, G, Bin; numerical update of
(eps_i, kappa) per coarse state. Labels S/T/W by fly-averaged movement fraction.
Rationale: everything the routing question needs is shared (precision from all flies); everything V2/V6
perturbed (rates, texture) has per-fly parameters to absorb fly differences instead of a fake state.

## Model comparison design (for the later validation, fixed now)
Held-out flies cannot be predicted (their E_i, eps_i are unknown); a held-out day for all flies has unknown
kappa. Use STAGGERED leave-one-day-out: fly i holds out day ((i mod D)); kappa for every day is then
estimated from the other flies, and each fly's parameters from its other days.

## Unit test (this step): 8 heterogeneous flies x 4 days from P3 with a day trend
Report: EM steps non-decreasing within 0.05 nats, fitted kappa_S vs truth, Xn_c/Yn_c vs truth, run time.
Informal; not a gate.
