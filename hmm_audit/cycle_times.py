# Descriptive (2026-11-03, Sheyum's question; no test): completed loops on each fly's minute symbols by loop erasure.
# A 3-state loop erased from the path is a completed cycle; forward = rotation of W->Q->S->W. Loop duration = minutes
# from entering the loop's first state to returning to it. Assigned to the phase of the completion minute.
import json, numpy as np, sqw
from sqw15 import seqs
FWD = {(1, 0, 2), (0, 2, 1), (2, 1, 0)}
def loops(y, ph):
    out = []; stack = []                      # entries [state, time entered]
    prev = None
    for t, x in enumerate(y):
        if x < 0: stack = []; prev = None; continue
        if x == prev: continue
        prev = x
        idx = next((k for k in range(len(stack)) if stack[k][0] == x), None)
        if idx is None: stack.append([x, t]); continue
        loop = [s for s, _ in stack[idx:]]
        if len(loop) == 3: out.append((ph[t], 1 if tuple(loop) in FWD else -1, t - stack[idx][1]))
        del stack[idx + 1:]; stack[idx][1] = t
    return out
O = {(d['g'], d['fly']): d for d in json.load(open(sqw.OUTD + 's_observable.json'))}
res = {}
for g in ['CS', 'CSDD', 'per0LD', 'per0DD', 'fumin', 'sss']:
    nP = 1 if g in ('CSDD', 'per0DD') else 2
    for p in range(nP):
        fp, fm, durf, durb, net, hrs = [], [], [], [], [], []
        for i in range(len(sqw.files(g))):
            if not O[(g, i)]['main'][p]['ok']: continue
            y, ph, _ = seqs(g, i); y = np.asarray(y); ph = np.asarray(ph)
            h = np.sum((ph == p) & (y >= 0)) / 60; L = [l for l in loops(y, ph) if l[0] == p]
            a = sum(1 for l in L if l[1] == 1); b = sum(1 for l in L if l[1] == -1)
            fp.append(a / h); fm.append(b / h); net.append((a - b) / h); hrs.append(h)
            durf += [l[2] for l in L if l[1] == 1]; durb += [l[2] for l in L if l[1] == -1]
        J = np.array([O[(g, i)]['main'][p]['J_chain'] or 0 for i in range(len(sqw.files(g))) if O[(g, i)]['main'][p]['ok']])
        res[f'{g}|{p}'] = dict(n=len(fp), fwd_per_h=float(np.mean(fp)), back_per_h=float(np.mean(fm)), net_per_h=float(np.mean(net)), J_chain_mean=float(J.mean()),
                              h_per_net_cycle=float(1 / np.mean(net)) if np.mean(net) > 0 else None,
                              dur_fwd_median_min=float(np.median(durf)) if durf else None, dur_back_median_min=float(np.median(durb)) if durb else None,
                              dur_fwd_mean_min=float(np.mean(durf)) if durf else None, dur_back_mean_min=float(np.mean(durb)) if durb else None,
                              n_loops_fwd=len(durf), n_loops_back=len(durb))
json.dump(res, open(sqw.OUTD + 'cycle_times.json', 'w'))
print('block     flies | completed loops per h: fwd  back  net (chain J) | h per NET cycle | loop duration fwd median/mean | back median/mean (min) | n loops f/b')
for k, v in res.items():
    hp = f"{v['h_per_net_cycle']:.1f}" if v['h_per_net_cycle'] else '  —'
    print(f"{k:9s} {v['n']:3d} | {v['fwd_per_h']:.3f} {v['back_per_h']:.3f} {v['net_per_h']:+.3f} ({v['J_chain_mean']:+.3f}) | {hp:>5s} | {v['dur_fwd_median_min']:.0f}/{v['dur_fwd_mean_min']:.0f} | {v['dur_back_median_min']:.0f}/{v['dur_back_mean_min']:.0f} | {v['n_loops_fwd']}/{v['n_loops_back']}")
