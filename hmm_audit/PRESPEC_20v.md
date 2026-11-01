# PRESPEC 20v — merged paper, test 2: waiting-time irreversibility beyond the observed chain
Date 2026-11-01. Written before anything runs; the go-ahead to RUN is Sheyum's (not given by this file).
Owner: Plan B feasibility session. Design agreed with the S/Q/W audit session. Results go to RESULTS_LEDGER_planB.md.

## Question
Do the dwell times in Q and W depend on where the fly came from and where it goes next, more than a first-order chain
allows? In a first-order Markov chain the dwell in a state is independent of its entry and exit neighbours, so any
dependence is irreversibility the transition counts miss. It is the waiting-time route to hidden entropy production
(Martinez, Bisker, Horowitz & Parrondo 2019 Nat Commun 10:3542; Skinner & Dunkel 2021 PRL 127:198101;
van der Meer, Ertel & Seifert 2022 PRX 12:031025; Harunari, Dutta, Polettini & Roldan 2022 PRX 12:041026).

## Data, symbols, blocks
1-min symbols exactly as sqw_ref.sqw_symbols (theta 0.5). DD: light_dark=False (one whole-record reference). LD:
light_dark=True. Fly-phases without a Q/W reference are excluded (20p amendment B).
PRIMARY: per0 DD whole record; CS DD whole record if uploaded before the run (Holm across the primary blocks present).
DESCRIPTIVE: CS, per0LD, fumin, sss, light and dark phases, with the 6-h-rate null (below), because within-phase drift
(20o) makes LD records non-stationary.

## Statistic (per fly x phase)
Runs of Q (and separately of W) with a known entry symbol e and exit symbol x (neighbouring minutes non-missing, in
the same phase). Only the mixed pairs carry direction:
- Q dwells: (e, x) = (W, S) vs its reverse (S, W).
- W dwells: (Q, S) vs (S, Q).
S dwells are NOT used (5-min floor); S enters only as an entry or exit event.
Dwell bins fixed now (minutes): 1, 2, 3, 4-5, 6-9, >= 10; pseudocount 0.5 per bin.
sigma_WT = sum over s in {Q, W} and the two mixed pairs of (n_{s,e,x} / hours) x KL( psi_{s,e,x} || psi_{s,x,e} ),
in nats/h, where psi is the binned dwell distribution and hours = minutes in the phase / 60.

## Null and excess (as 20p)
Null A: 200 records per fly-phase simulated from that fly's own mc_fit chain (12-h phases in LD; whole record in DD),
same length, missing minutes reinserted at their observed positions, scored identically.
LD descriptive null: the same from the 6-h-rate chain (mc_fit with phases(T, '6h')).
E_WT = sigma_WT(obs) - median sigma_WT(null); per-fly p = (1 + #{null >= obs}) / 201.

## Decision rule, fixed now
Per primary block: SUPPORTED if Holm p < 0.05 for the group Wilcoxon of E_WT against 0 with median E_WT > 0, AND the
fraction of flies with p < 0.05 exceeds Binomial(n, 0.05) at p < 0.05. Otherwise NOT SUPPORTED.
Wording if supported: "dwell times in quiet / active wake depend on the direction of passage beyond what the observed
chain allows: a waiting-time lower bound on hidden irreversibility of E_WT nats/h." E_WT is an EXCESS over the null,
not an EPR. The bound statement uses sigma_MC + E_WT, as 20p did with the triplet excess.
Wording if not: "no waiting-time evidence of hidden irreversibility at 1-min resolution." The coarse resolution is
stated (Q runs median 1-2 min), so a null result is weak.
sss gets one descriptive line: E_WT where the observed current is ~0 (C2) bears on a "hidden cycle without a
visible current" (Martinez et al. 2019). It is not part of the decision rule.

## Matched null, narrowest block, time box
Same fly, same length, same mask, same estimator. Narrowest block: one fly's whole DD record (LD: one fly's light or
dark minutes pooled over days). Expected run time under 1 h; no extension.
