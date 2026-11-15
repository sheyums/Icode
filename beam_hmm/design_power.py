"""Planning computation for the repeat experiment (DESIGN_NOTE.md). NOT an analysis of data: it uses only the
per-fly summary numbers already in RESULTS_LEDGER.md (pilot medians and fly-to-fly SDs) and simulates the
pre-specified one-sample Wilcoxon signed-rank test at the Holm worst-case level alpha/m.
Model: per-fly values ~ Normal(d, 1) with standardised effect d = pilot median / pilot SD. Two effect sizes per
measure: the pilot value, and half of it (pilot effects from ~10 flies, and those selected for being significant,
are optimistic: the 'winner's curse', e.g. Button et al. 2013 Nat Rev Neurosci 14:365).
Usage: python -I design_power.py <out.json>"""
import json, sys
import numpy as np
from scipy.stats import wilcoxon

MEASURES = [  # name, pilot median, pilot SD, Holm family size m, source
    ('A1 light D (end-resting asleep vs pause)', 0.293, 0.147, 2, 'ledger A1 weeks 2-3'),
    ('A4 walking speed dark - light', 0.231, 0.20, 2, 'ledger A2-A4 weeks 2-3'),
    ('A5 M1 sleep slope per week', 28.0, 24.75, 4, 'results/a5/a5.json'),
    ('A5 M3 light-D slope per week', 0.0310, 0.0609, 4, 'results/a5/a5.json'),
    ('A5 M4 walking-speed slope per week', 0.0726, 0.200, 4, 'results/a5/a5.json'),
    ('A5 M2 activity slope per week', 9.0, 33.6, 4, 'results/a5/a5.json'),
    ('A2 PRE - MID pause fraction (illustrative 0.03)', 0.03, 0.070, 2, 'SD from ledger A2; effect chosen, not observed'),
]
NS = [6, 8, 10, 12, 15, 20, 25, 30, 40, 50, 60, 80, 100]
REPS = 4000


def power(d, n, alpha, rng):
    x = rng.normal(d, 1.0, size=(REPS, n))
    p = np.array([wilcoxon(r, alternative='two-sided', method='exact' if n <= 50 else 'approx').pvalue for r in x])
    return float((p < alpha).mean())


def main(out):
    rng = np.random.default_rng(np.random.SeedSequence(20261115))
    res = []
    for name, med, sd, m, src in MEASURES:
        d = med / sd; row = dict(measure=name, pilot_median=med, pilot_sd=sd, holm_m=m, source=src, d_pilot=d)
        for tag, dd in (('pilot', d), ('half', d / 2)):
            curve = {}; n80 = None
            for n in NS:
                pw = power(dd, n, 0.05 / m, rng); curve[n] = pw
                if pw >= 0.8 and n80 is None:
                    n80 = n
                if n80 is not None and pw >= 0.95:
                    break
            row[tag] = dict(d=dd, power=curve, n_for_80pct=n80)
        res.append(row)
        print(f"{name:52s} d={d:5.2f}  n80 pilot={row['pilot']['n_for_80pct']}  n80 half={row['half']['n_for_80pct']}", flush=True)
    json.dump(dict(reps=REPS, alpha_family=0.05, ns=NS, rows=res), open(out, 'w'), indent=1)


if __name__ == '__main__':
    main(sys.argv[1])
