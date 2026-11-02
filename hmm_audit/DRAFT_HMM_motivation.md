# Draft text: motivating the hidden Markov model (2026-11-02; for the PNAS manuscript)

## Introduction (one sentence, where S/Q/W is introduced)
Because the three states are inferred from beam crossings with a threshold, and because sleep and wake may themselves
contain internal states that behaviour does not reveal directly, we also asked whether a hidden-state description
changes any conclusion.

## Results, first paragraph of "Is there anything beneath the observed cycle?"
We fitted each fly with a three-state hidden Markov model for three reasons. First, the symbols are scored, not
observed: a minute near the quiet/active threshold can be misclassified, and a model with an explicit emission matrix
asks whether the circulation survives once misclassification is allowed for. Second, coarse-graining can only lower
entropy production (Esposito, 2012), so dissipation measured on the observed states is a lower bound, and hidden states
could carry dissipation the observed cycle does not show. Third, fly sleep has internal structure, including deeper
stages within long bouts (van Alphen et al., 2013; Tainton-Heap et al., 2021) and covert processes visible only in
transition probabilities (Wiggin et al., 2020); a hidden-state model is the natural way to ask whether behaviour alone
reveals such states. The model's hidden states coincided with quiet wake, active wake and sleep (emission
probabilities 80-100% diagonal); its circulation agreed in sign with the observed one in 81-100% of flies, including
the reversal in sleepless by day; and it gave the same ranking of Onsager conductances (SI Table). Its gain in
likelihood over the observed chain was small and short-range: matched by one extra minute of memory, and BIC preferred
the observed chain in about three quarters of flies. We therefore report physical quantities from the observed chain,
whose entropy production is a rigorous lower bound, and use the hidden model as a check. Where the data do carry
information beyond the observed chain (sleep-length memory; direction-dependent dwell times) we measure it directly from
the data, without a model.

## Methods
Open "Hidden Markov model: robustness check and model-based estimate" with the three reasons above, one sentence each.

## Reason -> result map (for referees)
- misclassification -> sign agreement 81-100%; sleepless reversal holds in the HMM (fig_sss_reversal)
- hidden dissipation -> TUR (20q) and 3-min memory excess (20p), Fig A4 c-d
- internal states -> BIC vs chain (20n), memory gain = one extra minute (20o), sleep-length memory (C7, also at the
  hidden level), direction-dependent active-wake dwells (C9*)

## References
Esposito M (2012) Phys Rev E 85:041125. van Alphen B et al. (2013) J Neurosci 33:6917. Tainton-Heap LAL et al. (2021)
Curr Biol 31:578. Wiggin TD et al. (2020) PNAS 117:10024.

## SI table: Onsager L and effective affinity, observed vs hidden (computed 2026-11-02 from s_observable.json)
genotype | L observed [95% CI] | L hidden (counts on HMM paths) | L hidden (model rates) | A_eff obs | A_eff hidden | chain A
CS      | 0.272 [0.228,0.301] | 0.215 [0.172,0.228] | 0.235 [0.187,0.246] | 0.57 | 1.10 | 0.47
CS DD   | 0.262 [0.223,0.309] | (no HMM fits)       |                     | 0.51 |      | 0.50
per0 LD | 0.286 [0.255,0.316] | 0.245 [0.228,0.272] | 0.254 [0.229,0.291] | 0.43 | 0.95 | 0.42
per0 DD | 0.352 [0.294,0.399] | 0.330 [0.200,0.361] | 0.335 [0.215,0.368] | 0.64 | 0.97 | 0.64
fumin   | 0.243 [0.193,0.261] | 0.237 [0.171,0.276] | 0.218 [0.165,0.239] | 1.11 | 1.45 | 1.09
sss     | 0.145 [0.117,0.177] | 0.144 [0.119,0.178] | 0.113 [0.083,0.150] | 0.70 | 1.10 | 0.79
(L in cycles^2 nats^-1 h^-1, fly bootstrap 1000; A_eff = EPR/|J| median over fly-phases with J != 0, nats per cycle.)
