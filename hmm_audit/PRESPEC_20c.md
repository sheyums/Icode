# PRESPEC 20c — SECONDARY design: 6 h transition blocks (before running)
Date 2026-10-23. Requested by the user after the steps-per-fit extrapolation put 6 h blocks at routing
SD ~0.24-0.39 (borderline). Secondary to PRESPEC 20b; everything not stated here is as in 20b.

## Model
Per fly: shared emission matrix B; FOUR transition matrices, for transitions into minutes at ZT 0-6, 6-12,
12-18, 18-24. k = (K-1) + 4K(K-1) + 2K = 13, 32, 59 for K = 2, 3, 4. Genotypes: CS, per0LD, fumin (user's choice;
per0DD has no light phases). ALL flies of these genotypes are fitted, not a subset, so that no selection of
flies can follow results.

## Order (each step logged before the next)
1. Synthetic, 20 flies x 5760 min per case:
   - V1c: the V1b truth applied to 6 h blocks (blocks 1-2 = V1b light, 3-4 = V1b dark).
     Pass: K=3 chosen in >= 70%; |mean bias| < 0.1 and SD < 0.3 for Xn, Yn in every block.
   - V2c: the V2b truth. Pass: K=2 chosen in >= 80%.
2. Realistic truths (after the 20b primary fits exist): 5 non-degenerate flies per genotype (CS, per0LD,
   fumin), chosen with seed 21 from the 20b K=3 fits; each fly's fitted 20b model (its P_light for blocks 1-2,
   P_dark for 3-4, its B, its record length) is simulated 20 times and refitted with the 6 h model.
   Reported: SD and bias of Xn, Yn per block, per genotype. The 6 h design is called USABLE for a genotype if the
   median SD over its 5 flies is < 0.3 in every block.
3. Fly fits with the 6 h model for all flies of the three genotypes (theta = 0.5; K = 2, 3, 4; BIC as 20b),
   routing per block (Xn_pred, Yn_pred, FFBS counts attributed to the completion block), reported for every
   genotype, with the step-2 verdict printed beside it. A genotype judged not usable is still reported, marked
   as such.
