# Checks that sqw_ref.py reproduces the analysis scripts' outputs (20b, theta 0.5). Run from hmm_audit/.
import sys,os,json,numpy as np
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import sqw_ref as R, sqw, sqw2
from sqw2_common import em as em_pipe
fits={(d['g'],d['fly']):d for d in map(json.loads,open(sqw.OUTD+'fit2_0.5.jsonl'))}
RT=json.load(open(sqw.OUTD+'route2_0.5.json')); ok=True
def check(name,cond):
    global ok; ok&=bool(cond); print(('PASS ' if cond else 'FAIL ')+name,flush=True)
# 1. symbols, every fly
bad=0;n=0
for g in sqw.GROUPS:
    for i,p in enumerate(sqw.files(g)):
        x=R.read_counts(p); y,_=R.sqw_symbols(x,0.5,light_dark=(g!='per0DD')); y0,_,_=sqw2.seqs(g,i,0.5); n+=1; bad+=not np.array_equal(y,y0)
check(f'S/Q/W symbols identical for {n-bad}/{n} flies',bad==0)
# 2. em_step == batched EM (one start, 25 iterations, from the same start)
g,i=('CS',3); y,_=R.sqw_symbols(R.read_counts(sqw.files(g)[i]),0.5); ph=R.phases(len(y))
rng=np.random.default_rng(5); pi=rng.dirichlet(np.ones(3)); P=rng.dirichlet(np.ones(3),(2,3)); B=rng.dirichlet(np.ones(3),3)
p1,P1,B1=pi.copy(),P.copy(),B.copy()
for _ in range(25): p1,P1,B1,_=R.em_step(p1,P1,B1,y,ph)
import sqw2_common as C
Sn=1; pib=pi[None].copy(); Pb=P[None].copy(); Bb=B[None].copy()
for _ in range(25):
    E=C._E(Bb,y); A,c=C._forward(pib,Pb,E,ph); T=len(y); Bk=np.empty_like(A); Bk[-1]=1
    for t in range(T-2,-1,-1): Bk[t]=np.einsum('skl,sl->sk',Pb[:,ph[t+1]],E[t+1]*Bk[t+1])/c[t+1][:,None]
    G=A*Bk; G/=G.sum(2,keepdims=True); Rr=E*Bk/c[:,:,None]
    for p in range(2):
        idx=np.flatnonzero(ph[1:]==p)+1; Xi=np.einsum('tsk,tsl->skl',A[idx-1],Rr[idx])*Pb[:,p]; Pb[:,p]=Xi/Xi.sum(2,keepdims=True)
    oh=np.stack([(y==v) for v in range(3)],1).astype(float); Bn=np.einsum('tsk,tv->skv',G,oh); Bb=Bn/Bn.sum(2,keepdims=True); pib=G[0]
check(f'em_step x25 == batched EM update (max diff {max(abs(P1-Pb[0]).max(),abs(B1-Bb[0]).max()):.1e})',max(abs(P1-Pb[0]).max(),abs(B1-Bb[0]).max())<1e-10)
# 3. em_fit reproduces stored fits (same seed) for 2 flies, K=3
for g,i in (('CS',3),('sss',7)):
    y,_=R.sqw_symbols(R.read_counts(sqw.files(g)[i]),0.5); ph=R.phases(len(y)); seed=20+1000*sqw.GROUPS.index(g)+i
    ll=R.em_fit(y,ph,2,3,20,seed)[0]; check(f'em_fit {g} fly {i} logL {ll:.6f} == stored {fits[(g,i)]["F"]["3"]["ll"]:.6f}',abs(ll-fits[(g,i)]['F']['3']['ll'])<1e-6)
# 4. FFBS / Viterbi counts reproduce route2 (same seed)
for g,i in (('CS',3),('fumin',10),('per0DD',2)):
    d=fits[(g,i)]; y,_=R.sqw_symbols(R.read_counts(sqw.files(g)[i]),0.5,light_dark=(g!='per0DD'))
    ph=R.phases(len(y),'none' if g=='per0DD' else '12h'); nph=ph.max()+1; pi,P,B=(np.array(d[k]) for k in ('pi3','P3','B3'))
    Z=R.ffbs(pi,P,B,y,ph,200,20+1000*sqw.GROUPS.index(g)+i); hours=[np.sum(ph==p)/60 for p in range(nph)]
    S=R.posterior_summaries(Z,ph,nph,hours); rt=RT[f'{g}|{i}']
    check(f'{g} {i}: FFBS route counts == pipeline',np.allclose(S['route_counts'],rt['ffbs']))
    check(f'{g} {i}: FFBS pair counts == pipeline',np.allclose(S['pair_counts'],rt['pairs_ffbs']))
    zv=R.viterbi(pi,P,B,y,ph); check(f'{g} {i}: Viterbi route counts == pipeline',np.allclose(R.route_counts(zv,ph,nph),rt['vit']))
    # c_pred sign == J_model sign == affinity sign
    for p in range(nph):
        e=R.model_energetics(P[p]); c=R.c_pred(P[p])
        if abs(e['J'])<1e-12:   # a cycle edge is numerically zero: no cycle, sign undefined (PRESPEC 20e amendment A)
            check(f'{g} {i} phase {p}: no cycle (|J| {abs(e["J"]):.1e}, A {e["A"]:+.1f}) -- sign test not applicable',True); continue
        check(f'{g} {i} phase {p}: sign c_pred {c:+.3f} == sign J_model {e["J"]:+.4f} == sign A {e["A"]:+.3f}',np.sign(c)==np.sign(e['J'])==np.sign(e['A']))
print('ALL PASS' if ok else 'SOME CHECKS FAILED')
