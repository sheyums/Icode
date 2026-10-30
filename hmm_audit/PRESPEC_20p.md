# PRESPEC 20p — does the short-range memory carry irreversibility? Block KL divergence (before running)
Date 2026-10-30. Follows 20n/20o: the S/Q/W symbols are not first-order Markov, and the departure is short-range
(one extra minute of memory captures 1.25-2x the HMM's gain). Question: is that memory TIME-ASYMMETRIC, i.e. does the
observed sequence carry irreversibility beyond what its pair statistics (the observed first-order chain) show? If yes,
that is a model-free lower bound on hidden dissipation, and the physics motivation for a hidden layer
(Roldan & Parrondo 2010 PRL 105:150607; 2012 PRE 85:031129; Martinez et al. 2019 Nat Commun 10:3542).
Genotypes: all five (CS, fumin, per0LD, per0DD, sss) plus CS DD (new data, 25 flies, 4 days; user upload 2026-10-30;
checked DD: no activity step at the old light switches, -0.014 +- 0.008 vs CS LD +0.078). theta 0.5. Symbols exactly as
20b; CS DD handled exactly as per0DD (whole-record Q/W reference, one transition matrix).

## Estimator (per fly x phase; phase = light / dark for LD genotypes, whole record for DD)
- Blocks of n = 1..4 consecutive minutes, all in the same phase, none missing. Counts c_n(b), b in {Q,W,S}^n.
- p_n(b) = (c_n(b) + 0.5) / (N_n + 0.5 * 3^n)  (pseudocount on all 3^n blocks; impossible blocks and their reverses
  get equal mass and contribute 0).
- D_n = sum_b p_n(b) ln[p_n(b) / p_n(reverse b)];  D_1 = 0;  d_n = D_n - D_{n-1}, in nats/h (x60).
- d_2 is the observed first-order chain's irreversibility. For ANY stationary first-order chain d_3 = d_4 = d_2.
- PRIMARY statistic: Delta3 = d_3 - d_2.  SECONDARY: Delta4 = d_4 - d_3; and Delta3_QW, the same Delta3 computed only
  from blocks containing no S (wake-only), where the 5-min sleep criterion cannot act.

## Null (sample-size matched: every replicate goes through the same estimator on the same mask and phases)
- Null A (primary): 200 records simulated from the fly's fitted first-order chain (MC1, 20n procedure; 2 matrices LD,
  1 DD). Per fly-phase p_A = (1 + #{Delta3* >= Delta3_obs}) / 201 (one-sided: excess irreversibility).
- Null B (guard against within-phase drift, 20o): 200 records from MC1 with four 6-h matrices (LD and DD alike).
- Excess E = Delta3_obs - median(Delta3*_A), per fly-phase.
- Model predictions (not nulls): 50 records each from the fly's HMM (20b fit; CS DD fitted now with the 20b procedure,
  K = 3, 20 starts, seed 20 + 1000*5 + fly) and from MC2 (20o procedure). Report the median Delta3 each model predicts
  and the fraction of fly-phases whose Delta3_obs lies inside the model's central 95%.

## Tests and reading fixed now
- Per genotype x phase block (10 blocks): Wilcoxon signed-rank of E against 0 (two-sided), Holm across the 10 blocks;
  and the fraction of fly-phases with p_A < 0.05 against Binomial(n, 0.05).
- The memory carries EXTRA IRREVERSIBILITY in a block if Holm p < 0.05 with median E > 0, AND the fraction with p_A < 0.05
  exceeds the binomial (p < 0.05), AND the same holds under null B. Report d_3/d_2 there as the factor by which the
  observed chain underestimates irreversibility (a lower bound on that factor; d_4 also reported).
- The HMM ACCOUNTS for it if the HMM's predicted Delta3 brackets the observed in >= 80% of fly-phases of that block.
- If no block qualifies: the short memory is time-symmetric at the 3-4 min scale; this test gives no thermodynamic
  motivation for a hidden layer. That is a result and will be reported as such.
- Expected under the alternative "misclassification only": not determined in advance (noise on an irreversible hidden
  chain can produce some Delta3 > 0); the HMM-prediction comparison is what speaks to it.
- Limitations stated now: (i) a time-symmetric non-Markov feature (the >= 5-min sleep run criterion) could shift Delta3
  without hidden dissipation; Delta3_QW is the guard. (ii) Plug-in KL estimators are biased upward with block length;
  the matched null absorbs this, which is why no Delta3 is interpreted without it. (iii) KL irreversibility of
  behavioural symbols is a statistical arrow of time; calling it dissipation assumes the symbols are a coarse-graining of
  a Markov process whose states are even under time reversal.
- Narrowest block: one fly's light or dark minutes pooled over its days (DD: whole record).

## Amendment A (2026-10-30, written while the main run was in progress; NO result had been read)
User information: CS DD still oscillates (subjective day vs night), unlike per0DD, so a single whole-record phase is
non-stationary for CS DD. Added, for BOTH DD genotypes (per0DD as the matched comparison):
- The estimator run separately in the two 12-h clock halves (ZT 0-12, 12-24 of the recording clock; with a free-running
  period near 24 h these approximate subjective day/night to within ~1-2 h over 4 days; tau is not estimated).
- Null A2: 200 records from MC1 with two 12-h matrices (same halves). Null B unchanged (6-h matrices).
- Same tests and qualifying rule as the main plan (A2 in place of A), Holm across these 4 extra blocks.
- The whole-record DD rows of the main plan are still reported; for CS DD the halves rows are the ones to read.
No HMM / MC2 predictions for the halves (not needed for the qualifying rule).
