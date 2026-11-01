"""PRESPEC 20t F4, sleep lost per fly: window ZT15-24 (last 9 h of the night), deprivation day = day index 3
(located in planB_f3f4.log: group sleep 2-12 min/h in that window vs 38-50 on days 0-2); baseline = mean of days 0-2.
PASS: n >= 15 and P90/P10 of sleep lost >= 2."""
import json, os, numpy as np
H = np.array(json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'planB_f3f4.json')))['F4_hourly'])
w = lambda d: H[:, d*24+15:d*24+24].sum(1)
base = np.mean([w(0), w(1), w(2)], 0); dep = w(3); lost = base - dep
p10, p50, p90 = np.percentile(lost, [10, 50, 90])
print(f"F4 n {len(lost)} | baseline ZT15-24 sleep median {np.median(base):.0f} min | deprivation-night {np.median(dep):.0f} min "
      f"(range {dep.min():.0f}-{dep.max():.0f}) | sleep lost median {p50:.0f} min, P10 {p10:.0f}, P90 {p90:.0f}, "
      f"range {lost.min():.0f}-{lost.max():.0f} | P90/P10 {p90/p10:.2f} | CV {lost.std(ddof=1)/lost.mean():.2f}")
print('F4 PASS' if len(lost) >= 15 and p10 > 0 and p90 / p10 >= 2 else 'F4 FAIL')
print('per fly lost (min):', ' '.join(f'{v:.0f}' for v in np.sort(lost)))
