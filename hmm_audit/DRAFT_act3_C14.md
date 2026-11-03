# Draft Results text, act 3: sleep amount and the cycle's direction (C14)

Status: draft for the paper, 2026-11-03. The observed-level analysis is pre-registered (PRESPEC 22, 482898d). The
hidden-level panel is a post hoc HMM check and is labelled as one. Numbers are from RESULTS_LEDGER.md ("PRESPEC 22",
"POST HOC to PRESPEC 22", and the sensitivity note).

## Paragraph

Could the reversal in *sleepless* simply reflect how little these flies sleep? *sleepless* and *fumin* are the two
shortest sleepers in our panel, yet they circulate in opposite directions. To ask whether sleep amount sets the
direction, we related each fly's cycle to its own sleep within each genotype and phase (Fig. Xa,b). The direction per
transition did not depend on sleep amount in any genotype or phase (Spearman rho between -0.26 and +0.23, all
corrected p = 1; pooled within-genotype rho = 0.00, p = 0.98). Where the cycle flux did depend on sleep, it did so
through the number of transitions: *fumin* flies that slept more completed more forward cycles per hour at night
(rho = +0.59, corrected p = 0.003), with no change in direction per transition. At matched sleep in the dark (about
3 min h^-1 in both genotypes), *fumin* circulated forward and *sleepless* did not.

The hidden Markov model sharpens this (Fig. Xc,d). The direction per sleep onset was again independent of sleep in
every block. The flux scaled with sleep, and it did so with the sign of each genotype's own direction: more sleep meant
more forward cycles per hour in *fumin* at night (rho = +0.54), and more BACKWARD cycles per hour in *sleepless* by day
(rho = -0.51; the same with the fly's observed sleep, rho = -0.49). If sleeping little merely weakened an otherwise
forward cycle, the *sleepless* flies that sleep most should come closest to forward circulation. Instead they run the
cycle backward most often. Sleep amount sets how often the cycle turns; the genotype sets which way.

## Figure legend (draft)

Fig. X. Sleep amount sets how often the sleep-wake cycle turns, not its direction. Each point is one fly in one phase.
Open circles: light phase; filled circles: dark phase (or constant darkness); crosses: no net cycle. (a) Cycle flux J
of the observed chain, positive for forward circulation (active wake -> quiet wake -> sleep -> active wake).
(b) Direction per transition c_obs. (c, d) The same quantities from the hidden Markov model: flux averaged over 200
posterior paths, and net forward turns per sleep onset; x is the model's sleep occupancy. Panels a and b are
pre-registered; c and d are a post hoc check.

## Wording rules
- "post hoc" or "a check with the hidden model" for c and d; never "confirmed" or "pre-registered".
- Never "sleep amount does not affect the cycle": it affects the RATE (how often).
- The counterfactual sentence rests on sss light at the hidden level (corrected p 0.023 within the post hoc family);
  at the observed level the same sign is not significant (rho -0.26, p 0.11). Say so in the SI.
