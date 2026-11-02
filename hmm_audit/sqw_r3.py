# Recompute model-based quantities with the rule of three (2026-11-02): observed chain (mc_fit matrices per fly x
# phase, counts N from s_observable.json) and hidden model (20b K=3 matrices, mean posterior-path counts from
# route2_0.5.json). Writes sqw_out/r3_model.json and prints before/after group medians (infinities ranked highest
# for 'before', as reported in 20s-A).
import json, numpy as np, sqw, sqw2, sqw_ref as R
from sqw15 import seqs, DD
from sqw13 import mc_fit
O = json.load(open(sqw.OUTD + 's_observable.json')); F = sqw2.load(0.5); RT = json.load(open(sqw.OUTD + 'route2_0.5.json'))
out = []
for d in O:
    g, i = d['g'], d['fly']; nP = 1 if g in DD else 2
    y, ph, _ = seqs(g, i); _, _, Pm = mc_fit(y, ph, nP)
    hid = F.get((g, i)) if g != 'CSDD' else None
    rt = RT.get(f'{g}|{i}') if g != 'CSDD' else None
    for p in range(nP):
        r = d['main'][p]
        if not r['ok']: continue
        N = np.array(r['N'], float)
        e = R.model_energetics_r3(Pm[p], N); ls = R.landscape(R.rule_of_three_matrix(Pm[p], N)[0])
        row = dict(g=g, fly=i, p=p, cycle=(r['J_chain'] is not None and abs(r['J_chain']) >= 1e-12),
                   A_old=r['A_chain'], EPR_old=r['EPR_chain'], A_r3=float(e['A']), EPR_r3=float(e['EPR']), r3=bool(e['r3']),
                   E21_r3=float(ls['E2_minus_E1']), E31_r3=float(ls['E3_minus_E1']), F_SW_r3=float(ls['F_SW']), F_WS_r3=float(ls['F_WS']))
        if hid is not None and not hid['deg']:
            Ph = np.array(hid['P3'])[p]; Nh = np.mean(np.array(rt['pairs_path'])[:, p], 0)
            eh = R.model_energetics_r3(Ph, Nh)
            row.update(EPR_hmm_model_old=r.get('EPR_hmm_model'), EPR_hmm_model_r3=float(eh['EPR']), A_hmm_r3=float(eh['A']), r3_hmm=bool(eh['r3']))
        out.append(row)
json.dump(out, open(sqw.OUTD + 'r3_model.json', 'w'))
def med(v):
    v = np.array([np.inf if (x is not None and not np.isfinite(x) and not np.isnan(x)) else x for x in v if x is not None], float)
    v = v[~np.isnan(v)]; return np.median(v) if len(v) else np.nan
print('genotype phase | n cycle | r3 used | A median old -> r3 | chain EPR median old -> r3 | HMM model EPR median old -> r3 (r3 used)')
for g in ['CS', 'CSDD', 'per0LD', 'per0DD', 'fumin', 'sss']:
    for p in ([0] if g in DD else [0, 1]):
        rr = [x for x in out if x['g'] == g and x['p'] == p and x['cycle']]
        hh = [x for x in out if x['g'] == g and x['p'] == p and 'EPR_hmm_model_r3' in x]
        print(f"{g:7s} {p} | {len(rr):2d} | {sum(x['r3'] for x in rr):2d} | {med([x['A_old'] for x in rr]):+.2f} -> {med([x['A_r3'] for x in rr]):+.2f} | "
              f"{med([x['EPR_old'] for x in rr]):.3f} -> {med([x['EPR_r3'] for x in rr]):.3f} | "
              + (f"{med([x['EPR_hmm_model_old'] for x in hh]):.3f} -> {med([x['EPR_hmm_model_r3'] for x in hh]):.3f} ({sum(x['r3_hmm'] for x in hh)})" if hh else '—'))
nonfin = [x for x in out if x['cycle'] and not (np.isfinite(x['A_r3']) and np.isfinite(x['EPR_r3']))]
print('observed fly-phases with a cycle still non-finite after r3:', len(nonfin))
