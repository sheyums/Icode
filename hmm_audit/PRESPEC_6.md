# PRESPEC 6 — how much of the multi-phase bout structure is pooling heterogeneity? (before running)
Date 2026-10-19. Prompted by PRESPEC 5 P1 (sleep K >= 2 in every genotype). No analysis below run yet.

## Question
Is the pooled K=2-vs-K=1 gain larger than what between-WINDOW (hence between-fly and
time-of-day) heterogeneity of single-exponential bouts would produce?

## Data
- Bouts exactly as PRESPEC 5 (interior runs, sleep xmin = 300 s, wake xmin = 2 s, dt = 1 s).
- Block = fly x window (the narrowest block; also where each HMM was fitted). A bout belongs
  to the window it starts in. Windows with < 3 bouts of the class are excluded (count reported);
  the pooled observed statistic is recomputed on the SAME retained bouts.
- Primary: sleep. Secondary: wake. All 4 genotypes.

## Null H0 (heterogeneity only)
Within each window, excess e = d - (xmin - 1) >= 1 is geometric with a WINDOW-SPECIFIC
p_w = 1/mean(e_w) (MLE), n_w bouts as observed. Plug-in rates from ~5-10 bouts are noisy, so
H0 carries MORE heterogeneity than the truth: the null is conservative toward "it is all
heterogeneity". Heterogeneity WITHIN a window (drift inside 4 h) is not tested by this.

## Statistics (each computed identically on observed and on every replicate)
- T1 = pooled dL = logL(K=2) - logL(K=1), discrete left-truncated hyperexponential
  (pooled.py fit(), 10 starts for K=2).
- T2 = within-window overdispersion V = sum_w (n_w - 1)[s_w^2/m_w^2 - (1 - 1/m_w)] / sum_w (n_w - 1),
  e-values, s_w^2 unbiased. 0 in expectation-ish for geometric; > 0 = heavier than exponential
  INSIDE a window. Its small-n bias is handled by the matched simulation, not analytically.
- 500 replicates, seed default_rng(6). Report observed, null mean and 95% range, one-sided
  p = (1 + #null >= obs)/(1 + 500), and the share explained = null-mean T1 / observed T1.
## Read-out fixed now
- p(T1) < 0.05 and p(T2) < 0.05: multi-phase structure exceeds between-window heterogeneity
  and is present within windows.
- p(T1) >= 0.05: the pooled K >= 2 is not distinguishable from between-window heterogeneity.
