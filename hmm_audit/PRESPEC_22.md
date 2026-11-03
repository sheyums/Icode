# PRESPEC 22 — does a fly's cycle flux depend on how much it sleeps?

Date 2026-11-03. Question from Sheyum. This is an exception to the analysis freeze, which Sheyum asked for. Written and
committed BEFORE any per-fly sleep amount has been set against any per-fly cycle quantity. Already known, and so NOT
evidence: the GROUP-level facts. sss and fumin are the shortest sleepers. sss reverses the cycle by day. fumin runs it
forward. A pooled across-genotype correlation is therefore largely fixed in advance by the group means, and is
DESCRIPTIVE only. The test is WITHIN genotype.

Numbering: 22 is used here; the Sh^mns and extra-silencing PRESPECs proposed on 2026-11-02 become 23 and 24.

## Data and quantities (no new definitions)
- Genotypes as Sheyum named them: CS (LD), per0 DD, fumin, sss. CS DD and per0 LD are reported as SECONDARY blocks,
  descriptive only, and are not in the Holm family.
- Blocks are fly x phase (12-h light/dark). per0 DD is one whole-record block. That gives 7 primary blocks: CS L,
  CS D, per0 DD, fumin L, fumin D, sss L, sss D. No re-binning, no ZT sub-windows.
- Sleep = minutes scored S per hour of non-missing minutes in the block (sqw15.seqs symbols, theta 0.5, as 20s-A).
- J = observed-chain cycle flux J_chain (cycles/h, s_observable.json, fitted matrix; NOT changed by the rule of
  three), with the sign convention of fig_sss_reversal: positive = forward W->Q->S->W.
- c_obs = the minute-level current of 20d (P(last wake minute before onset is Q) - P(first wake minute after waking
  is Q)). It is a direction PER TRANSITION, so it does not scale with how many transitions a fly makes.
- Only blocks with ok references are used (as 20s-A).

## Why two outcomes (a mechanical link fixed in advance)
|J| cannot exceed the rate of sleep onsets, so a fly that sleeps little MUST have a small |J|. A correlation of sleep
with |J| or with J is partly built in. The question that is not built in is the DIRECTION. So:
- PRIMARY: Spearman rho(sleep, c_obs) within each of the 7 blocks.
- CO-PRIMARY: Spearman rho(sleep, J) within each block, over fly-phases with a cycle (|J_chain| >= 1e-12), as in
  the sss reversal. A positive rho here is reported together with the c_obs result, and is never read alone as
  "sleep sets the direction".
- Null: a within-block permutation of sleep against the outcome, 10000 permutations, seed 22. This null has the
  block's own n by construction (rule 5). Two-sided p. Holm over the 14 tests (7 blocks x 2 outcomes).

## Secondary, descriptive
- No-cycle fly-phases (|J_chain| < 1e-12): their median sleep against the median sleep of fly-phases with a cycle,
  per block. Are the no-cycle flies the shortest sleepers?
- Pooled figure: J (and c_obs) against sleep, every fly-phase, coloured by genotype, light open and dark filled,
  with per-block Spearman lines. Descriptive: genotype is a confounder.
- Stratified summary: the n-weighted mean of the within-block Fisher z of rho(sleep, c_obs) over the 7 blocks, with
  a permutation p from permuting within blocks. It answers "does sleep predict direction beyond genotype".
- Overlap check: are fumin and sss sleep ranges overlapping? Where they do, does genotype or sleep separate c_obs?
  Descriptive.

## Decision rules (fixed now)
- "Sleep amount predicts the cycle's direction within a genotype": rho(sleep, c_obs) Holm < 0.05 in that block.
- The act 3 claim ("sleep amount and the breaking of detailed balance are dissociable") is WITHDRAWN if both of
  these hold:
  (a) the stratified Fisher-z summary for c_obs has p < 0.05 with a positive sign (more sleep, more forward); and
  (b) in the fumin/sss overlap range, sleep and not genotype orders c_obs.
  If (a) holds without (b), the claim is REWORDED to "sleep amount predicts direction within genotypes, but does not
  explain the fumin/sss contrast".
- A J-only correlation, with c_obs null, is reported as the mechanical rate link, not as direction.
- No other genotype, phase split, sleep measure (e.g. bout number or bout length) or outcome is added after
  looking. Any of them would be a separate, labelled post hoc.

## Output
`sqw22.py` -> `sqw_out/s22_sleep_vs_J.json` and `figures/fig22_J_vs_sleep.png`. Ledger entry with commit hashes.
