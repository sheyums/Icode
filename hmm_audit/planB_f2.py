"""PRESPEC 20t F2: separability of clock and homeostat in real flies.
r = multiple correlation of H(t) (Process S from the fly's own S minutes) on {cos wt, sin wt}, tau = 24 h.
Primary (tau_r, tau_d) = (6 h, 3 h); range over the 36-point grid reported. Pass per genotype: r < 0.8 in > half of flies.
Also a data-identity check: CS c_obs medians vs 20s-A (light +0.042, dark +0.147).
Usage: python3 -I planB_f2.py <data_root>   (one subfolder per genotype)"""
import sys, os, glob, json, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import sqw_ref as R
from planB_f1f5 import h_path, multiple_r, GRID

SETS = {'CS': ('CS_activity_files_1sec_res', True), 'per0LD': ('per0LD_activity_1sec_res', True),
        'fumin': ('Fmn_activity_files_1sec', True), 'sss': ('sleepless_activity_1sec_res', True),
        'per0DD': ('per0DD_activity_files_1sec_res', False)}

def files(root, sub):
    return sorted(f for f in glob.glob(os.path.join(root, sub, '**', '*'), recursive=True) if os.path.isfile(f))

def main(root):
    out = {}
    for g, (sub, ld) in SETS.items():
        rows = []
        for f in files(root, sub):
            y, refs = R.sqw_symbols(R.read_counts(f), theta=0.5, light_dark=ld)
            r0 = multiple_r(h_path(y, 6, 3))
            rg = [multiple_r(h_path(y, a, b)) for a in GRID for b in GRID]
            rec = dict(file=os.path.basename(f), r=r0, rmin=min(rg), rmax=max(rg), fracS=float(np.mean(y == R.S)))
            if g == 'CS':
                ph = R.phases(len(y)); rec['c_obs'] = R.c_obs_minutes(y, ph, 2).tolist()
            rows.append(rec)
        r = np.array([x['r'] for x in rows])
        frac = float(np.mean(r < 0.8))
        out[g] = dict(n=len(rows), median_r=float(np.median(r)), frac_below_08=frac, PASS=frac > 0.5,
                      rmax_grid_median=float(np.median([x['rmax'] for x in rows])),
                      frac_below_08_at_grid_max=float(np.mean([x['rmax'] < 0.8 for x in rows])), flies=rows)
        print(f"{g:7s} n {len(rows)} | r(6h,3h) median {np.median(r):.2f} IQR [{np.percentile(r,25):.2f},{np.percentile(r,75):.2f}] "
              f"max {r.max():.2f} | r<0.8 in {frac:.2f} -> {'PASS' if frac > 0.5 else 'FAIL'} | grid: median of per-fly max r "
              f"{out[g]['rmax_grid_median']:.2f}, r<0.8 at worst grid point in {out[g]['frac_below_08_at_grid_max']:.2f}")
        if g == 'CS':
            c = np.array([x['c_obs'] for x in rows])
            print(f"   identity check CS c_obs median light {np.nanmedian(c[:,0]):+.3f} dark {np.nanmedian(c[:,1]):+.3f} (20s-A: +0.042, +0.147)")
    json.dump(out, open(os.path.join(HERE, 'planB_f2.json'), 'w'))

if __name__ == '__main__':
    main(sys.argv[1])
