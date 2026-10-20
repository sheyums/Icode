# PRESPEC 10, C3 only, resumable: identical to counts_analysis.C3job (same code, same seeds);
# each finished (genotype, fly, phase) is appended to C3_results.jsonl and skipped on restart.
import os; os.environ['OMP_NUM_THREADS']='1'
import sys, json; sys.path.insert(0,'/tmp/claude-0/-home-user-Icode/1f385af9-cc34-5eb3-921b-4648ebd049a8/scratchpad/hmm_audit')
from counts_analysis import *
OUT='C3_results.jsonl'
def done_keys():
    if not os.path.exists(OUT): return set()
    return {(r['g'],r['fly'],r['phase']) for r in map(json.loads,open(OUT))}
if __name__=='__main__':
    data={g:load(g) for g in GEN}; dk=done_keys()
    jobs=[(g,fi,p,x) for g in GEN for fi,x in enumerate(data[g]) for p in [0,1] if (g,fi+1,p) not in dk]
    print('remaining jobs',len(jobs),flush=True)
    with Pool(4) as pl:
        for r in pl.imap_unordered(C3job,jobs):
            with open(OUT,'a') as f: f.write(json.dumps(r,default=float)+'\n')
    R3=[json.loads(l) for l in open(OUT)]
    for g in GEN:
        for p in [0,1]:
            r=[x for x in R3 if x['g']==g and x['phase']==p]
            deg=np.mean([min(x['L1_12'],x['L1_13'])<0.05 for x in r]); unid=np.mean([(x['Xn_range']>0.2) or (x['Yn_range']>0.2) for x in r])
            print(f"{g:6s} {'light' if p==0 else 'dark '} flies {len(r)} | H1 degenerate (min L1<0.05) {deg:.2f} | median L1(H1,H2) {np.median([x['L1_12'] for x in r]):.3f} L1(H1,H3) {np.median([x['L1_13'] for x in r]):.3f} | route-unidentified (restart range>0.2) {unid:.2f} | median restarts within 0.5 nats {np.median([x['n_near'] for x in r]):.0f}")
            print(f"        median Xn {np.nanmedian([x['Xn'] for x in r]):+.3f} Yn {np.nanmedian([x['Yn'] for x in r]):+.3f} | median Xn_range {np.nanmedian([x['Xn_range'] for x in r]):.3f} Yn_range {np.nanmedian([x['Yn_range'] for x in r]):.3f}")
            ho=np.array([[x['heldout_ll_per_s'][str(K)] if str(K) in x['heldout_ll_per_s'] else x['heldout_ll_per_s'][K] for K in [2,3,4]] for x in r if 'heldout_ll_per_s' in x])
            if len(ho): print(f"        held-out log-lik per s, median over flies: K2 {np.median(ho[:,0]):.5f} K3 {np.median(ho[:,1]):.5f} K4 {np.median(ho[:,2]):.5f}; flies where K3>K2 {np.mean(ho[:,1]>ho[:,0]):.2f}, K4>K3 {np.mean(ho[:,2]>ho[:,1]):.2f}",flush=True)
