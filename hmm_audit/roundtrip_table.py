# Per-fly round-trip table (PRESPEC 20g), built only from sqw_ref.py functions and the stored 20b fits.
# Output: roundtrip_per_fly_theta<theta>.csv, one row per fly x phase x level. Usage: python3 roundtrip_table.py [theta]
import sys,os,json,csv,numpy as np
from multiprocessing import Pool
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import sqw_ref as R, sqw, sqw2
TH=float(sys.argv[1]) if len(sys.argv)>1 else 0.5
FITS={(d['g'],d['fly']):d for d in map(json.loads,open(sqw.OUTD+f'fit2_{TH}.jsonl'))}
def one(a):
    g,i=a; path=sqw.files(g)[i]; x=R.read_counts(path); ld=g!='per0DD'; nph=2 if ld else 1
    y,_=R.sqw_symbols(x,TH,light_dark=ld); ph=R.phases(len(y),'12h' if ld else 'none'); hours=[np.sum(ph==p)/60 for p in range(nph)]
    rows=[]; Co=R.round_trips_observed(x,TH,ld)
    levels=[('observed',[Co[p] for p in range(nph)])]
    d=FITS.get((g,i))
    if d is not None and 'P3' in d and not d['deg']:
        pi,P,B=(np.array(d[k]) for k in ('pi3','P3','B3')); Z=R.ffbs(pi,P,B,y,ph,200,7000+1000*sqw.GROUPS.index(g)+i)
        Ha=np.stack([R.round_trips_hidden(z,ph,nph) for z in Z]); Hm=np.stack([R.round_trips_hidden(z,ph,nph,5,10) for z in Z])
        levels+= [('hidden_all',[Ha[:,p] for p in range(nph)]),('hidden_matched',[Hm[:,p] for p in range(nph)])]
    for lev,cs in levels:
        for p in range(nph):
            s=R.round_trip_stats(cs[p],hours[p])
            if s is None: continue
            rows.append([g,i,os.path.basename(path),'whole' if not ld else ['light','dark'][p],lev,TH,round(s['n'],2)]
                        +[round(v,4) for v in s['frac']]+[round(v,4) for v in s['rate']]+[round(s['p_q_entry'],4),round(s['p_q_exit'],4)]
                        +[round(v,4) for v in s['dependence']]+[round(s['frac_paths_used'],3)])
    return rows
if __name__=='__main__':
    with Pool(4) as pl: out=pl.map(one,[(g,i) for g in sqw.GROUPS for i in range(len(sqw.files(g)))],chunksize=2)
    hdr=['genotype','fly_index','file','phase','level','theta','n_round_trips']+[f'frac_{t}' for t in R.ROUND_TRIPS]+[f'rate_per_h_{t}' for t in R.ROUND_TRIPS]\
        +['p_quiet_entry','p_quiet_exit']+[f'dependence_{t}' for t in R.ROUND_TRIPS]+['frac_paths_used']
    fn=f'roundtrip_per_fly_theta{TH}.csv'
    with open(fn,'w',newline='') as f:
        w=csv.writer(f); w.writerow(hdr); [w.writerow(r) for rows in out for r in rows]
    print('wrote',fn,sum(len(r) for r in out),'rows')
