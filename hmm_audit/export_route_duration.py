# Export per-bout route/duration records (PRESPEC 20k) for the user's MATLAB JSD analysis.
# Observed: one row per bridged sleep bout with neighbours on both sides. Hidden: one row per H3 run on each of the
# first 20 posterior paths per fly (path index given). Durations in minutes; routes: 1 = quiet (Q / H1),
# 0 = active (A / H2), empty = not classified (neighbouring wake period too short).
import json,csv,sys; sys.path.insert(0,'.'); import sqw
res=json.load(open('sqw_out/k_route_duration.json'))
files={(g,i):sqw.files(g)[i].split('/')[-1] for g in sqw.GROUPS for i in range(len(sqw.files(g)))}
def b(v): return '' if v is None else int(bool(v))
hdr=['genotype','fly_index','file','phase','sleep_min','wake_before_min','wake_after_min','entry_quiet','exit_quiet']
with open('route_duration_observed_theta0.5.csv','w',newline='') as f:
    w=csv.writer(f); w.writerow(hdr)
    for r in res:
        for x in r['obs']:
            w.writerow([r['g'],r['fly'],files[(r['g'],r['fly'])],'whole' if r['g']=='per0DD' else ['light','dark'][x['phase_on']],
                        round(x['sleep_min'],3),round(x['wake_before_min'],3),round(x['wake_after_min'],3),b(x['q_in']),b(x['q_out'])])
with open('route_duration_hidden_theta0.5.csv','w',newline='') as f:
    w=csv.writer(f); w.writerow(hdr[:3]+['path']+hdr[3:])
    for r in res:
        for k,path in enumerate(r['hidrec']):
            for x in path:
                w.writerow([r['g'],r['fly'],files[(r['g'],r['fly'])],k,'whole' if r['g']=='per0DD' else ['light','dark'][x['phase_on']],
                            x['sleep_min'],x['wake_before_min'],x['wake_after_min'],b(x['q_in']),b(x['q_out'])])
print('exported')
