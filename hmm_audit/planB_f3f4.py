"""PRESPEC 20t F3 (clock-precision bound informativeness) and F4 (deprivation data: format, n, spread).
Usage: python3 -I planB_f3f4.py <CS_DD dir> <per0DD dir> <deprivation dir>
F3: per 24-h window the circular-mean (centroid) clock time of activity counts; intervals between successive markers;
    per-fly mean interval removed; pooled s^2 = SS / (n_int - n_fly); fly bootstrap (2000) of s; bound B = 2 tau / s^2.
    Rhythmic: chi-square periodogram (Sokolove & Bushell 1978) on 1-min activity counts, periods 20-28 h in 6-min steps,
    peak Qp above chi2 0.99 quantile with P-1 df. PASS: >= 15 of 25 CS DD flies rhythmic AND bound CI upper/lower < 4.
    s from RHYTHMIC flies (primary; the PRESPEC counts rhythmic flies) and from all flies (sensitivity).
F4: sleep per fly per day per hour (sqw_ref.sleep_bouts, 5-min criterion); the deprivation window is located from the
    data (user: last 9 h of night 3, i.e. ZT15-24 of day index 2); sleep lost = mean baseline-day sleep in the window
    minus deprivation-day sleep in the window. PASS: readable, n >= 15, P90/P10 of sleep lost >= 2."""
import sys, os, glob, json, numpy as np
from scipy.stats import chi2
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import sqw_ref as R

def files(d): return sorted(f for f in glob.glob(os.path.join(d, '**', '*'), recursive=True) if os.path.isfile(f))

def periodogram(m):
    """m: 1-min counts. Returns (best period min, Qp, threshold)."""
    best = (np.nan, -np.inf, np.nan); X = m - m.mean(); den = np.sum(X**2)
    for P in range(1200, 1681, 6):
        K = len(m) // P; x = m[:K * P].reshape(K, P)
        Mh = x.mean(0); N = K * P
        Qp = K * N * np.sum((Mh - x.mean())**2) / np.sum((x - x.mean())**2)
        thr = chi2.ppf(0.99, P - 1)
        if not np.isfinite(best[2]) or Qp - thr > best[1] - best[2]:   # peak = largest margin above threshold
            best = (P, Qp, thr)
    return best

def markers(m):
    nd = len(m) // 1440; out = []
    for d in range(nd):
        w = m[d * 1440:(d + 1) * 1440]; th = 2 * np.pi * np.arange(1440) / 1440
        z = np.sum(w * np.exp(1j * th)); out.append(np.angle(z) % (2 * np.pi) * 1440 / (2 * np.pi))
    return np.array(out)

def intervals(mk):
    d = np.diff(mk); d = (d + 720) % 1440 - 720   # wrapped difference in (-12, 12] h
    return 1440 + d                                # minutes

def pooled_s(iv_list):
    ss = sum(np.sum((iv - iv.mean())**2) for iv in iv_list); df = sum(len(iv) for iv in iv_list) - len(iv_list)
    tau = np.mean(np.concatenate(iv_list))
    return np.sqrt(ss / df) if df > 0 else np.nan, tau

def f3(d, label, rng):
    rows = []
    for f in files(d):
        x = R.read_counts(f); m = x[:len(x) // 60 * 60].reshape(-1, 60).sum(1).astype(float)
        P, Qp, thr = periodogram(m); iv = intervals(markers(m))
        rows.append(dict(file=os.path.basename(f), period_h=P / 60, Qp=Qp, thr=thr, rhythmic=bool(Qp > thr), iv=iv.tolist()))
    res = dict(n=len(rows), n_rhythmic=int(sum(r['rhythmic'] for r in rows)))
    for key, sel in (('rhythmic', [r for r in rows if r['rhythmic']]), ('all', rows)):
        ivs = [np.array(r['iv']) for r in sel]
        if len(ivs) < 2: res[key] = None; continue
        s, tau = pooled_s(ivs)
        bs = []
        for b in range(2000):
            idx = rng.integers(len(ivs), size=len(ivs)); sb, tb = pooled_s([ivs[i] for i in idx]); bs.append((sb, tb))
        bs = np.array(bs); B = 2 * (tau / 60) / (s / 60)**2
        Bb = 2 * (bs[:, 1] / 60) / (bs[:, 0] / 60)**2
        lo, hi = np.percentile(Bb, [2.5, 97.5]); slo, shi = np.percentile(bs[:, 0] / 60, [2.5, 97.5])
        res[key] = dict(n_flies=len(ivs), n_int=int(sum(len(v) for v in ivs)), tau_h=tau / 60, s_h=s / 60,
                        s_CI=[slo, shi], B_per_h=B, B_CI=[lo, hi], CI_ratio=hi / lo)
        print(f"F3 {label:7s} [{key:8s}] flies {len(ivs)} intervals {res[key]['n_int']} | tau {tau/60:.2f} h | s {s/60:.2f} h "
              f"[{slo:.2f}, {shi:.2f}] | B = 2 tau/s^2 = {B:.1f} k_B/h [{lo:.1f}, {hi:.1f}] | CI ratio {hi/lo:.2f}")
    print(f"F3 {label:7s} rhythmic {res['n_rhythmic']} of {res['n']} | periods (rhythmic) median "
          f"{np.median([r['period_h'] for r in rows if r['rhythmic']]) if res['n_rhythmic'] else float('nan'):.2f} h")
    res['flies'] = rows
    return res

def f4(d):
    rows = []; hourly = []
    for f in files(d):
        x = R.read_counts(f); T = len(x); sl = np.zeros(T, bool)
        for s, n in R.sleep_bouts(x, True): sl[s:s + n] = True
        nh = T // 3600; hs = sl[:nh * 3600].reshape(nh, 3600).sum(1) / 60.0   # sleep min per hour
        hourly.append(hs); rows.append(dict(file=os.path.basename(f), days=T / 86400))
    nh = min(len(h) for h in hourly); H = np.array([h[:nh] for h in hourly])
    print(f"F4 files {len(rows)} | length {rows[0]['days']:.4f} d | hours {nh}")
    print("F4 group mean sleep (min per hour) by day (rows) x ZT hour (cols):")
    for dd in range(nh // 24):
        print(f"   day {dd}: " + ' '.join(f"{v:4.0f}" for v in H[:, dd * 24:(dd + 1) * 24].mean(0)))
    if nh % 24: print(f"   day {nh // 24} (partial): " + ' '.join(f"{v:4.0f}" for v in H[:, (nh // 24) * 24:].mean(0)))
    return H, rows

def main(csdd, per0dd, depr):
    rng = np.random.default_rng(2026); out = {}
    out['F3_CSDD'] = f3(csdd, 'CS DD', rng); out['F3_per0DD'] = f3(per0dd, 'per0 DD', rng)
    r = out['F3_CSDD']; k = r.get('rhythmic')
    ok = r['n_rhythmic'] >= 15 and k is not None and k['CI_ratio'] < 4
    print(f"F3 criterion: rhythmic {r['n_rhythmic']} >= 15 and CI ratio {k['CI_ratio'] if k else float('nan'):.2f} < 4 -> {'PASS' if ok else 'FAIL'}")
    H, rows = f4(depr)
    out['F4_hourly'] = H.tolist(); out['F4_files'] = rows
    json.dump(out, open(os.path.join(HERE, 'planB_f3f4.json'), 'w'))

if __name__ == '__main__':
    main(*sys.argv[1:4])
