"""DESCRIPTIVE figure for PRESPEC 21 (C13): pooled 1-s activity profiles at the two ends of wake bouts.
m_start(tau): active fraction tau seconds after waking (start of the wake bout).
m_end(tau):   active fraction tau seconds before sleep onset (end of the wake bout), i.e. the end window time-reversed.
d(tau) = m_start - m_end per fly; median over flies with a fly-bootstrap 95% band (2000 replicates, seed 2026).
Profiles: mean over flies of each fly's mean over its bouts, fly-bootstrap 95% bands. Same bouts as planB_21
(wake bouts between consecutive sleep bouts, >= 2L; flies with < 5 bouts excluded). Both pre-declared settings.
Usage: python3 -I fig_planB_21.py <data root> <CS DD folder>  -> figures/fig_planB21_profiles_{primary,supp}.{png,pdf}"""
import sys, os, json, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import sqw_ref as R
from planB_21 import traces, SETTINGS
BLUE, ORANGE, INK, INK2, GRID, SURF = '#2a78d6', '#eb6834', '#0b0b0b', '#52514e', '#e4e3df', '#fcfcfb'
NAMES = {'per0DD': 'per0 DD', 'CSDD': 'CS DD', 'CS': 'CS LD', 'per0LD': 'per0 LD', 'fumin': 'fumin LD', 'sss': 'sss LD'}

def fly_profiles(x, L, b):
    a = (np.asarray(x) > 0).astype(float); sb = R.sleep_bouts(x, True); st = []; en = []
    for k in range(len(sb) - 1):
        w0 = sb[k][0] + sb[k][1]; w1 = sb[k + 1][0]
        if w1 - w0 < 2 * L: continue
        st.append(a[w0:w0 + L].reshape(-1, b).mean(1)); en.append(a[w1 - L:w1][::-1].reshape(-1, b).mean(1))
    return (np.mean(st, 0), np.mean(en, 0), len(st)) if len(st) >= 5 else None

def summarise(P, rng):
    S = np.array([p[0] for p in P]); E = np.array([p[1] for p in P]); D = S - E; n = len(P)
    bs = [rng.integers(n, size=n) for _ in range(2000)]
    band = lambda f: np.percentile([f(i) for i in bs], [2.5, 97.5], axis=0)
    return dict(n=n, ms=S.mean(0), me=E.mean(0), ms_b=band(lambda i: S[i].mean(0)), me_b=band(lambda i: E[i].mean(0)),
                d=np.median(D, 0), d_b=band(lambda i: np.median(D[i], 0)))

def style(ax):
    ax.set_facecolor(SURF)
    for s in ('top', 'right'): ax.spines[s].set_visible(False)
    for s in ('left', 'bottom'): ax.spines[s].set_color(INK2)
    ax.tick_params(colors=INK2, labelsize=8); ax.grid(True, color=GRID, lw=0.6); ax.set_axisbelow(True)

def panel_row(axs, s300, s120, name, first):
    for (ax_p, ax_d), s, (L, b) in ((axs[0:2], s300, SETTINGS['L300']), (axs[2:4], s120, SETTINGS['L120'])):
        t = (np.arange(L // b) + 0.5) * b
        for m, bd, col, lab in ((s['ms'], s['ms_b'], BLUE, 'after waking'), (s['me'], s['me_b'], ORANGE, 'before sleep onset (reversed)')):
            ax_p.fill_between(t, bd[0], bd[1], color=col, alpha=0.18, lw=0); ax_p.plot(t, m, color=col, lw=2, label=lab)
        lo_, hi_ = ax_p.get_ylim(); gap = 0.07 * (hi_ - lo_); ya, yb = s['ms'][-1], s['me'][-1]
        if abs(ya - yb) < gap:                                   # keep the two end labels apart (no collision)
            mid = (ya + yb) / 2; sgn = 1 if ya >= yb else -1; ya, yb = mid + sgn * gap / 2, mid - sgn * gap / 2
        ax_p.text(t[-1], ya, ' after waking', color=INK2, fontsize=7, va='center', ha='left')
        ax_p.text(t[-1], yb, ' before onset', color=INK2, fontsize=7, va='center', ha='left')
        ax_d.fill_between(t, s['d_b'][0], s['d_b'][1], color=BLUE, alpha=0.18, lw=0); ax_d.plot(t, s['d'], color=BLUE, lw=2)
        ax_d.axhline(0, color=INK2, lw=0.8)
        for ax in (ax_p, ax_d): style(ax); ax.set_xlim(0, L * 1.18 if ax is ax_p else L)
        ax_p.set_ylabel(f"{name} (n {s['n']})\nactive fraction", fontsize=8, color=INK)
        ax_d.set_ylabel('d = after − before', fontsize=8, color=INK)
        if first:
            ax_p.set_title(f'profiles, L {L} s, {b}-s bins', fontsize=9, color=INK)
            ax_d.set_title(f'median d(τ), L {L} s', fontsize=9, color=INK)

def make(groups, S, fname, title):
    fig, axes = plt.subplots(len(groups), 4, figsize=(15, 2.9 * len(groups) + 0.8), squeeze=False)
    fig.patch.set_facecolor(SURF)
    for r, g in enumerate(groups):
        panel_row(axes[r], S['L300'][g], S['L120'][g], NAMES[g], r == 0)
    for ax in axes[-1]: ax.set_xlabel('τ, seconds from the sleep boundary', fontsize=8, color=INK)
    h, l = axes[0][0].get_legend_handles_labels()
    fig.legend(h, l, loc='upper right', fontsize=8, frameon=False, ncol=2)
    fig.suptitle(title, fontsize=10, color=INK, x=0.01, ha='left')
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    os.makedirs(os.path.join(HERE, 'figures'), exist_ok=True)
    for ext in ('png', 'pdf'): fig.savefig(os.path.join(HERE, 'figures', f'{fname}.{ext}'), dpi=200, facecolor=SURF)
    plt.close(fig)

def main(root, csdd):
    T = traces(root, csdd); S = {}; rng = np.random.default_rng(2026)
    for st, (L, b) in SETTINGS.items():
        S[st] = {}
        for g, X in T.items():
            P = [p for p in (fly_profiles(x, L, b) for x in X) if p is not None]; S[st][g] = summarise(P, rng)
    make(['per0DD', 'CSDD'], S, 'fig_planB21_profiles_primary',
         'DESCRIPTIVE (PRESPEC 21 / C13): 1-s activity at the two ends of wake bouts, DD. Bands: fly bootstrap 95%.')
    make(['CS', 'per0LD', 'fumin', 'sss'], S, 'fig_planB21_profiles_supp',
         'DESCRIPTIVE, supplementary (PRESPEC 21 / C13): LD genotypes. Bands: fly bootstrap 95%.')
    out = {st: {g: {k: (v.tolist() if isinstance(v, np.ndarray) else v) for k, v in S[st][g].items()} for g in S[st]} for st in S}
    json.dump(out, open(os.path.join(HERE, 'fig_planB21_profiles.json'), 'w'))
    for st in S:
        for g in ('per0DD', 'CSDD'):
            s = S[st][g]
            print(f"{st} {g:7s} n {s['n']} | m_start first/last bin {s['ms'][0]:.3f}/{s['ms'][-1]:.3f} | m_end first/last bin {s['me'][0]:.3f}/{s['me'][-1]:.3f}")

if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
