# Round-trip composition per genotype x phase: observed (1-s counts, PRESPEC 20g) and hidden (posterior paths,
# definition matched to observed: amendment B, post hoc). theta 0.5. Means over flies.
import json,numpy as np,matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.family':'STIXGeneral','mathtext.fontset':'stix','font.size':8,'axes.linewidth':0.6})
INK='#1a1a1a'
S=json.load(open('sqw_out/g_summary.json')); M=json.load(open('sqw_out/g_matched.json'))
BLOCKS=[('CS','light'),('CS','dark'),('per0LD','light'),('per0LD','dark'),('per0DD','whole'),('fumin','light'),('fumin','dark'),('sss','light'),('sss','dark')]
NAME={'CS':'CS','per0LD':'$per^0$ LD','per0DD':'$per^0$ DD','fumin':'$fmn$','sss':'$sss$'}
COLS=['#2a78d6','#eb6834','#b9b8b3','#e4e3df']; LAB=['forward  Q→S→A  (H1→H3→H2)','reverse  A→S→Q  (H2→H3→H1)','direct  A→S→A','quiet both  Q→S→Q']
def hid(g,bn):
    p=0 if bn in ('light','whole') else 1; rows=[v[p] for k,v in M.items() if k.startswith(g+'|') and v[p] is not None]
    return np.array([r['frac'] for r in rows]).mean(0)
fig,axs=plt.subplots(1,2,figsize=(7.0,3.1),sharey=True)
ylab=[f"{NAME[g]}{'' if bn=='whole' else (' L' if bn=='light' else ' D')}" for g,bn in BLOCKS]; yy=np.arange(len(BLOCKS))[::-1]
for ax,(title,get) in zip(axs,[('(a) observed behaviour',lambda g,bn:np.array(S[f'0.5|observed|{g}|{bn}']['mean'])),('(b) hidden states (matched definition)',hid)]):
    for y,(g,bn) in zip(yy,BLOCKS):
        f=get(g,bn); left=0
        for j in range(4):
            ax.barh(y,f[j],left=left,height=0.72,color=COLS[j],edgecolor='white',lw=0.8); 
            if f[j]>=0.08: ax.text(left+f[j]/2,y,f'{f[j]*100:.0f}',ha='center',va='center',fontsize=6,color='white' if j<2 else INK)
            left+=f[j]
    ax.set_xlim(0,1); ax.set_xticks([0,0.25,0.5,0.75,1]); ax.set_xticklabels(['0','25','50','75','100%'])
    ax.set_title(title,fontsize=8,loc='left'); ax.tick_params(axis='y',length=0)
    for s in ('top','right'): ax.spines[s].set_visible(False)
axs[0].set_yticks(yy); axs[0].set_yticklabels(ylab)
axs[0].set_xlabel('share of sleep bouts'); axs[1].set_xlabel('share of sleep bouts')
h=[plt.Rectangle((0,0),1,1,color=c) for c in COLS]
fig.legend(h,LAB,loc='lower center',ncol=4,frameon=False,fontsize=6.5,handlelength=1.2,columnspacing=1.2)
fig.tight_layout(rect=(0,0.08,1,1),w_pad=1.0)
fig.savefig('figures/fig8_roundtrips.pdf'); fig.savefig('figures/fig8_roundtrips.png',dpi=300)
