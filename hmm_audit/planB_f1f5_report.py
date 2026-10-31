"""Report for PRESPEC 20t F1 / F5 from planB_f1f5.json. Criteria exactly as in PRESPEC_20t.md."""
import json, os, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
D = json.load(open(os.path.join(HERE, 'planB_f1f5.json')))
fl = D['flies']; meta = D['meta']
C4, C2 = 9.49, 5.99
ROW = {0: 'Q', 1: 'W', 2: 'S'}
MOD = ['02', '21']     # Q->S, S->W
def by(s): return [f for f in fl if f['scen'] == s]
def det(fs, kind, i, crit): return np.mean([f['LR_' + kind][i] > crit for f in fs])
print(f"seconds {D['seconds']:.0f}")
for s, m in meta.items():
    print(f"{s:11s} days {m['days']} a {m['a']} Eb {m['Eb']} b {m['b']:.3f} sdH_calib {m['sdH_calib']:.4f} Hc {m['Hc']:.3f}")
print('\nDetection rates (LRT > nominal chi2 5%): clock rows Q/W/S | homeostat rows Q/W/S | N0-null 95th pct LR')
res = {}
for s in meta:
    fs = by(s)
    c = [det(fs, 'clock', i, C4) for i in range(3)]; h = [det(fs, 'homeo', i, C2) for i in range(3)]
    res[s] = (c, h)
    print(f"{s:11s} n {len(fs)} | clock {c[0]:.2f} {c[1]:.2f} {c[2]:.2f} | homeo {h[0]:.2f} {h[1]:.2f} {h[2]:.2f}")
for s in ('N0', 'N0-2d'):
    fs = by(s)
    print(f"  {s} 95th pct LR: clock " + ' '.join(f"{np.percentile([f['LR_clock'][i] for f in fs], 95):.2f}" for i in range(3))
          + ' | homeo ' + ' '.join(f"{np.percentile([f['LR_homeo'][i] for f in fs], 95):.2f}" for i in range(3)))

print('\nEstimates per modulated transition (median over flies [IQR])')
est = {}
for s in meta:
    fs = by(s); m = meta[s]
    for tr in MOD:
        Ebt = np.array([f['coef_true'][tr][3] * f['sd_Htrue'] for f in fs])
        Ebh = np.array([f['coef'][tr][3] * f['sd_Hhat'] for f in fs])
        at = np.array([np.hypot(*f['coef_true'][tr][1:3]) for f in fs])
        ah = np.array([np.hypot(*f['coef'][tr][1:3]) for f in fs])
        pt = np.array([np.arctan2(f['coef_true'][tr][2], f['coef_true'][tr][1]) for f in fs])
        phh = np.array([np.arctan2(f['coef'][tr][2], f['coef'][tr][1]) for f in fs])
        perr = np.abs((phh - pt + np.pi) % (2 * np.pi) - np.pi) * 24 / (2 * np.pi)
        line = f"{s:11s} {tr}"
        if m['Eb'] > 0:
            ratio = Ebh / Ebt; sign = np.mean(np.sign(Ebh) == np.sign(Ebt))
            est[(s, tr)] = (np.median(ratio), sign)
            line += f" | Eb true {np.median(Ebt):+.3f} hat {np.median(Ebh):+.3f} ratio {np.median(ratio):.2f} [{np.percentile(ratio,25):.2f},{np.percentile(ratio,75):.2f}] sign ok {sign:.2f}"
        else:
            line += f" | Eb hat (truth 0) median {np.median(Ebh):+.3f} IQR [{np.percentile(Ebh,25):+.3f},{np.percentile(Ebh,75):+.3f}]"
        if m['a'] > 0:
            line += f" | a true {np.median(at):.2f} hat {np.median(ah):.2f} ratio {np.median(ah/at):.2f} | phase err median {np.median(perr):.1f} h"
        else:
            line += f" | a hat (truth 0) {np.median(ah):.2f}"
        print(line)
    print(f"{'':11s} tau_r true {np.median([f['tr_true'] for f in fs]):.1f} h hat median {np.median([f['tau_r'] for f in fs]):.0f} | "
          f"tau_d true {np.median([f['td_true'] for f in fs]):.1f} hat median {np.median([f['tau_d'] for f in fs]):.0f} | "
          f"r(H,clock) true {np.median([f['r_true'] for f in fs]):.2f} hat {np.median([f['r_hat'] for f in fs]):.2f} | frac S {np.median([f['frac_S'] for f in fs]):.2f}")

print('\nF5 thermodynamics (nats/h; median over flies)')
f5 = {}
for s in ('N0', 'DD-CH', 'DD-CH-half', 'LD2-CH'):
    fs = by(s)
    hk_t = np.array([f['hk_true'] for f in fs]); hk_h = np.array([f['hk_hat'] for f in fs])
    na_t = np.array([f['na_true'] for f in fs]); na_h = np.array([f['na_hat'] for f in fs])
    f5[s] = dict(hk=np.median(hk_h / hk_t), na=np.median(na_h / na_t) if s != 'N0' else np.nan,
                 na_t=np.median(na_t), na_h=np.median(na_h))
    print(f"{s:11s} hk true {np.median(hk_t):.4f} hat {np.median(hk_h):.4f} ratio {np.median(hk_h/hk_t):.2f} | "
          f"na true {np.median(na_t):.5f} hat {np.median(na_h):.5f} ratio {np.median(na_h/na_t):.2f} [{np.percentile(na_h/na_t,25):.2f},{np.percentile(na_h/na_t,75):.2f}] | "
          f"na share true {np.median(na_t/(na_t+hk_t)):.4f}")

print('\nF1 criteria')
ok_all = True
for s, n0 in (('DD-CH', 'N0'), ('LD2-CH', 'N0-2d')):
    c, h = res[s]
    k1 = min(c[0], c[2], h[0], h[2]) >= 0.80
    cn, hn = res[n0]
    k2 = c[1] <= 0.10 and h[1] <= 0.10 and max(cn) <= 0.10 and max(hn) <= 0.10
    k4 = all(0.75 <= est[(s, tr)][0] <= 1.25 and est[(s, tr)][1] >= 0.80 for tr in MOD)
    print(f"{s}: (1) detection min {min(c[0], c[2], h[0], h[2]):.2f} -> {k1} | (2) FP W-row clock {c[1]:.2f} homeo {h[1]:.2f}, "
          f"{n0} max clock {max(cn):.2f} homeo {max(hn):.2f} -> {k2} | (4) Eb ratio/sign "
          + ', '.join(f"{tr} {est[(s, tr)][0]:.2f}/{est[(s, tr)][1]:.2f}" for tr in MOD) + f" -> {k4}")
    ok_all &= k1 and k2 and k4
cC, hC = res['DD-C']; cH, hH = res['DD-H']
k3 = max(hC[0], hC[2]) <= 0.10 and max(cH[0], cH[2]) <= 0.10
print(f"(3) cross-talk: homeostat in DD-C Q/S {hC[0]:.2f}/{hC[2]:.2f}; clock in DD-H Q/S {cH[0]:.2f}/{cH[2]:.2f} -> {k3}")
ok_all &= k3
print('F1 PASS' if ok_all else 'F1 FAIL')
ok5 = True
for s in ('DD-CH', 'LD2-CH'):
    a = 0.80 <= f5[s]['hk'] <= 1.25; b = 0.5 <= f5[s]['na'] <= 2.0
    print(f"F5 {s}: hk ratio {f5[s]['hk']:.2f} -> {a} | na ratio {f5[s]['na']:.2f} -> {b}"); ok5 &= a and b
fl_ = f5['N0']['na_h'] / f5['DD-CH']['na_t']
print(f"F5 noise floor: N0 fitted na {f5['N0']['na_h']:.5f} / DD-CH true na {f5['DD-CH']['na_t']:.5f} = {fl_:.2f} -> {fl_ < 0.25}")
ok5 &= fl_ < 0.25
print('F5 PASS' if ok5 else 'F5 FAIL')
