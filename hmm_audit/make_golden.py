# Golden test vectors for the MATLAB port: fixed inputs (synthetic, no fly data) and sqw_ref.py outputs.
# All indices in the JSON are 0-BASED as in Python; PORTING_GUIDE.md says where to add 1.
import json,os,numpy as np
import sqw_ref as R
OUT=os.path.join(os.path.dirname(os.path.abspath(__file__)),'golden'); os.makedirs(OUT,exist_ok=True)
def clean(o):
    """JSON-safe: arrays -> lists; non-finite floats -> the strings 'Inf', '-Inf', 'NaN' (see PORTING_GUIDE)."""
    if isinstance(o,dict): return {k:clean(v) for k,v in o.items()}
    if isinstance(o,(list,tuple)): return [clean(v) for v in o]
    if isinstance(o,np.ndarray): return clean(o.tolist())
    if isinstance(o,(bool,np.bool_)): return bool(o)
    if isinstance(o,(int,np.integer)): return int(o)
    if isinstance(o,(float,np.floating)):
        o=float(o); return o if np.isfinite(o) else ('NaN' if np.isnan(o) else ('Inf' if o>0 else '-Inf'))
    return o
def dump(name,d): json.dump(clean(d),open(os.path.join(OUT,name),'w'),indent=0,allow_nan=False)
# G1: bridging edge cases (hand-built 1-s traces; 1 = active second)
cases={}
def trace(*segs):   # segs: (value, length)
    return np.concatenate([np.full(n,v) for v,n in segs]).astype(int)
cases['301_1_300_bridged']=trace((1,5),(0,301),(1,1),(0,300),(1,5))
cases['301_1_100_not_bridged']=trace((1,5),(0,301),(1,1),(0,100),(1,5))
cases['300_2_300_not_bridged']=trace((1,5),(0,300),(1,2),(0,300),(1,5))
cases['299_1_400_not_bridged']=trace((1,5),(0,299),(1,1),(0,400),(1,5))
cases['chain_300_1_300_1_300']=trace((1,5),(0,300),(1,1),(0,300),(1,1),(0,300),(1,5))
cases['at_start_no_left_flank']=trace((0,400),(1,1),(0,400),(1,5))
cases['exactly_300']=trace((1,3),(0,300),(1,3))
cases['299']=trace((1,3),(0,299),(1,3))
dump('G1_sleep_bouts.json',{k:dict(x=v,bouts_bridged=R.sleep_bouts(v,True),bouts_unbridged=R.sleep_bouts(v,False)) for k,v in cases.items()})
# G2: a synthetic 2-day count trace -> bouts, references, symbols (theta 0.5, 12h)
rng=np.random.default_rng(1); P0=np.array([[0.67,0.165,0.165],[0.05765,0.93,0.01235],[0.02074,0.01926,0.96]])
B0=np.array([[0.80,0.15,0.05],[0.15,0.85,0.00],[0.03,0.00,0.97]])
z,_=R.simulate(R.stationary(P0),np.stack([P0,P0]),B0,R.phases(2880),rng)
lam=np.array([0.05,0.6,0.0])          # per-second activity probability by hidden state (synthetic)
x=np.concatenate([(rng.random(60)<lam[k]).astype(int)*rng.integers(1,4,60) for k in z])
bouts=R.sleep_bouts(x,True); ev=R.wake_bout_events(x,bouts); y,refs=R.sqw_symbols(x,0.5)
dump('G2_symbols.json',dict(x=x,theta=0.5,bouts=bouts,n_wake_events=len(ev),refs=refs,y=y,
     note='x: 1-s counts from ZT0; y: 0=Q 1=W 2=S -1=missing; bouts: (start 0-based, length)'))
# G3: forward log-likelihood, one EM step, Viterbi, on fixed parameters and the G2 symbols
ph=R.phases(len(y)); pi=np.array([0.2,0.5,0.3])
P=np.stack([P0,np.array([[0.6,0.2,0.2],[0.05,0.9,0.05],[0.02,0.03,0.95]])]); B=B0*0.98+0.02/3
y3=y.copy(); y3[[5,6,700]]=-1
ll=R.loglik(pi,P,B,y3,ph); pi1,P1,B1,ll0=R.em_step(pi,P,B,y3,ph); zv=R.viterbi(pi,P,B,y3,ph)
dump('G3_hmm.json',dict(y=y3,ph=ph,pi=pi,P=P,B=B,loglik=ll,em_step=dict(pi=pi1,P=P1,B=B1,loglik_input=ll0),
     viterbi=zv,n_params={f'K{K}_nP{n}':R.n_params(K,n) for K in (2,3,4) for n in (1,2,4)},
     bic=dict(ll=-2000.0,K=3,nP=2,n=2880,value=R.bic(-2000.0,3,2,2880))))
# G4: labelling, degeneracy, routing, energetics on fixed matrices
Bs=np.array([[0.03,0.00,0.97],[0.80,0.15,0.05],[0.15,0.85,0.00]]); Ps=np.stack([P0[[2,0,1]][:,[2,0,1]]])
Pl,Bl,pil=R.label_states(Ps,Bs,np.array([0.5,0.2,0.3]))
rev=np.array([[0.67,0.165,0.165],[0.04125,0.93,0.02875],[0.033,0.023,0.944]])
mats={'V1_light':P0,'reversible':rev,'one_way_edge':np.array([[0.7,0.3,0.0],[0.1,0.85,0.05],[0.05,0.0,0.95]])}
dump('G4_matrices.json',dict(label=dict(P_in=Ps,B_in=Bs,P_out=Pl,B_out=Bl,pi_out=pil),
     degenerate=dict(B_distinct=B0,distinct=R.is_degenerate(B0),B_copy=np.array([[0.1,0.9,0],[0.12,0.88,0],[0,0,1.0]]),
                     copy=R.is_degenerate(np.array([[0.1,0.9,0],[0.12,0.88,0],[0,0,1.0]]))),
     matrices={k:dict(P=m,xn=R.xn_pred(m),yn=R.yn_pred(m),c=R.c_pred(m),pi=R.stationary(m),**R.model_energetics(m)) for k,m in mats.items()}))
# G5: counting on fixed hidden paths
zz=np.array([2,2,0,1,1,0,2,2,1,1,2,0,0,2,0,1,1,2,2,0,2,1]); phz=np.array([0]*11+[1]*11)
N=R.pair_counts(zz,phz,2); rc=R.route_counts(zz,phz,2)
Ncases={'all_positive':np.array([[0,12,30],[25,0,8],[15,20,0]]),'one_way_k5':np.array([[0,5,0],[0,0,4],[6,4,0]]),
        'one_way_k2':np.array([[0,2,0],[1,0,3],[4,3,0]]),'all_zero':np.zeros((3,3))}
dump('G5_counts.json',dict(z=zz,ph=phz,route_events=R.route_events(zz),route_counts=rc,pair_counts=N,
     xy_counts=R.xy_counts(np.array([[5,3,2,6],[1,1,4,0]])),
     epr_bound={k:dict(N=v,hours=10.0,epr=R.epr_count_bound(v,10.0)[0],substituted=R.epr_count_bound(v,10.0)[1],J=R.j_cycle_count(v,10.0)) for k,v in Ncases.items()},
     c_obs_minutes=dict(y=np.array([1,1,0,2,2,2,0,1,1,1,2,2,1,0,0,1,0,2]),ph=np.array([0]*9+[1]*9),
                        value=R.c_obs_minutes(np.array([1,1,0,2,2,2,0,1,1,1,2,2,1,0,0,1,0,2]),np.array([0]*9+[1]*9),2))))
# G6: FFBS posterior means with many paths (statistical check: MATLAB result must agree within the given SE)
Z=R.ffbs(pi,P,B,y3,ph,4000,7); S=R.posterior_summaries(Z,ph,2,[np.sum(ph==p)/60 for p in range(2)])
rp=np.stack([R.route_counts(z,ph,2) for z in Z])
dump('G6_ffbs.json',dict(note='inputs = G3; 4000 paths; compare MATLAB posterior mean route counts within 4 SE',
     mean_route_counts=rp.mean(0),se_route_counts=rp.std(0)/np.sqrt(len(Z)),c=S['c'],Xn=S['Xn'],Yn=S['Yn'],J=S['J'],EPR_bound=S['EPR_bound']))
print('golden written:',sorted(os.listdir(OUT)))
# G7: observable pair counts and round trips (PRESPEC 20f/20g)
g2=json.load(open(os.path.join(OUT,'G2_symbols.json'))); x2=np.array(g2['x']); y2=np.array(g2['y']); ph2=R.phases(len(y2))
zr=np.array([1]*12+[2]*6+[1]*3+[0]*2+[2]*8+[0]*11+[1]*4+[2]*7+[1]*15+[0]*1+[2]*9+[0]*12+[2]*3+[1]*20)
phr=np.zeros(len(zr),int); phr[len(zr)//2:]=1
Hp=np.stack([R.round_trips_hidden(zr,phr,2)[0],R.round_trips_hidden(zr,phr,2)[0][::-1]])
dump('G7_roundtrips.json',dict(
    observable_pair_counts=dict(note='inputs: G2 y, phases 12h',value=R.observable_pair_counts(y2,ph2,2)),
    round_trips_observed=dict(note='inputs: G2 x, theta 0.5, light_dark true',value=R.round_trips_observed(x2,0.5,True)),
    hidden_path=dict(z=zr,ph=phr,all=R.round_trips_hidden(zr,phr,2),matched_5_10=R.round_trips_hidden(zr,phr,2,5,10)),
    stats_single=dict(counts=[6,2,9,3],hours=12.0,**{k:v for k,v in R.round_trip_stats([6,2,9,3],12.0).items()}),
    stats_paths=dict(counts=Hp,hours=12.0,**{k:v for k,v in R.round_trip_stats(Hp,12.0).items()})))
print('G7 written')
# G8: dwell-time functions (PRESPEC 20h)
s8=np.array([2,2,2,1,1,0,2,2,-1,1,1,1,1,2,2,2,2]); p8=np.array([0]*12+[1]*5)
dump('G8_dwell.json',dict(seq=s8,ph=p8,edges=R.DWELL_EDGES,
    runs={nm:R.runs_censored(s8,p8,v,2) for nm,v in {'S':{2},'wake':{0,1},'Q':{0},'W':{1}}.items()},
    runs_no_phase_censoring=R.runs_censored(s8,p8,{0,1},2,censor_at_phase=False),
    km=dict(d=[3,2,1,4],e=[1,0,1,1],out=R.km_survival([3,2,1,4],[1,0,1,1])),
    km_weighted=dict(d=[3,2,1,4],e=[1,0,1,1],w=[0.5,1,1,2],out=R.km_survival([3,2,1,4],[1,0,1,1],[0.5,1,1,2])),
    hazard=dict(d=[3,2,1,4,7,15],e=[1,0,1,1,1,0],edges=[1,2,3,5,8,13,21],out=R.life_table_hazard([3,2,1,4,7,15],[1,0,1,1,1,0],[1,2,3,5,8,13,21])),
    geometric=dict(p=0.9,t=[0,1,5],out=R.geometric_survival(0.9,[0,1,5]))))
print('G8 written')
