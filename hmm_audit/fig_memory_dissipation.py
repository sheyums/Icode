# Explainer figure for PRESPEC 20p: does the short memory carry dissipation? Panels a-c are schematic; panel d is data.
import os,sys,json,numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Circle, FancyBboxPatch
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import sqw
COL={'Q':'#1baf7a','W':'#eb6834','S':'#2a78d6'}; INK='#0b0b0b'; INK2='#52514e'; GRID='#e4e3df'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.edgecolor':INK2,'axes.labelcolor':INK,'xtick.color':INK2,'ytick.color':INK2})
fig=plt.figure(figsize=(12,9.2),facecolor='white')
# ---------- (a) movie forwards and backwards ----------
ax=fig.add_axes([0.03,0.57,0.46,0.38]); ax.set_xlim(-2.6,20); ax.set_ylim(-0.6,4.2); ax.axis('off')
seq=list('WWWQQQSSSSSSWQWWWQQQ')
def strip(y0,s,label):
    for k,c in enumerate(s):
        ax.add_patch(FancyBboxPatch((k+0.06,y0),0.88,0.7,boxstyle='round,pad=0,rounding_size=0.08',fc=COL[c],ec='white',lw=1.5))
        ax.text(k+0.5,y0+0.35,c,ha='center',va='center',color='white',fontsize=8,fontweight='bold')
    ax.text(-0.2,y0+0.35,label,ha='right',va='center',color=INK,fontsize=9)
strip(2.9,seq,'forward'); strip(1.6,seq[::-1],'reversed')
ax.text(-2.6,3.95,'a   The same record played forwards and backwards (illustration)',fontsize=10.5,fontweight='bold',color=INK)
ax.add_patch(FancyBboxPatch((0.02,2.84),2.96,0.82,boxstyle='round,pad=0,rounding_size=0.1',fc='none',ec=INK,lw=1.4,ls='--'))
ax.text(1.5,2.55,'3-min window: W W W',ha='center',fontsize=7.5,color=INK2)
ax.text(0,0.95,'Pairs (2-min windows): every Q→W is eventually matched by a W→Q,\nso pair counts can look balanced, as if nothing runs in a preferred direction.',fontsize=8.5,color=INK,va='top')
ax.text(0,0.0,'Triplets (3-min windows) also see the ORDER inside a run: e.g. "W Q Q" (slowing down)\nvs "Q Q W" (abrupt start). If one is more common than its reverse, the movie has an\narrow of time that pairs cannot show. That extra arrow is a lower bound on dissipation.',fontsize=8.5,color=INK,va='top')
# ---------- (b) mechanism: hidden ring, lumped observation ----------
ax=fig.add_axes([0.52,0.57,0.46,0.38]); ax.set_xlim(-1.5,4.2); ax.set_ylim(-1.75,1.8); ax.set_aspect('equal'); ax.axis('off')
ax.text(-1.5,1.62,'b   Why it happens: a hidden cycle seen through a lumped camera',fontsize=10.5,fontweight='bold',color=INK)
pos={'a':(-0.9,0.35),'b':(0.0,-0.95),'c':(0.9,0.35)}; disp={'a':'Q','b':'Q','c':'W'}
for k,(x,y) in pos.items():
    ax.add_patch(Circle((x,y),0.36,fc=COL[disp[k]],ec='white',lw=2)); ax.text(x,y,k,ha='center',va='center',color='white',fontsize=12,fontweight='bold')
    ax.text(x,y-0.52 if k!='b' else y-0.5,f'shows as {disp[k]}',ha='center',va='top',fontsize=7.5,color=INK2)
for s,t in (('a','b'),('b','c'),('c','a')):
    (x0,y0),(x1,y1)=pos[s],pos[t]; ax.add_patch(FancyArrowPatch((x0,y0),(x1,y1),connectionstyle='arc3,rad=-0.28',arrowstyle='-|>',mutation_scale=16,lw=2,color=INK,shrinkA=17,shrinkB=17))
ax.text(0,1.05,'hidden circulation  a → b → c → a',ha='center',va='center',fontsize=8.5,color=INK2)
ax.text(1.55,-0.35,'The camera cannot tell a from b.\nObserved: … W Q Q W Q Q W …\n• Q→W and W→Q pairs balance\n• but inside each Q run the fly\n  always goes a then b, so runs\n  "end the same way" more than\n  they "start" that way:\n  triplets are asymmetric',fontsize=7.8,color=INK,va='center')
# ---------- (c) the measurement ----------
ax=fig.add_axes([0.08,0.17,0.36,0.30])
n=np.array([2,3,4]); null=np.array([0.10,0.15,0.20]); obs=np.array([0.10,0.20,0.25])
ax.plot(n,null,'-o',color=INK2,lw=2,ms=8,mec='white',mew=2,label='simulated from the observed first-order chain\n(no memory): rise from finite data only')
ax.plot(n,obs,'-o',color=COL['W'],lw=2,ms=8,mec='white',mew=2,label='real record')
ax.annotate('',xy=(3,0.198),xytext=(3,0.152),arrowprops=dict(arrowstyle='<->',color=INK,lw=1.3))
ax.text(2.35,0.215,'excess E\n= irreversibility the\nmemory adds',fontsize=8,color=INK,va='bottom',ha='left')
ax.set_xticks(n); ax.set_xticklabels(['pairs\n(2 min)','triplets\n(3 min)','quadruplets\n(4 min)'])
ax.set_ylabel('irreversibility seen with n-minute windows\n(k$_B$ per hour; schematic)'); ax.set_yticks([])
ax.set_xlim(1.7,4.4); ax.set_ylim(0.06,0.29); ax.grid(axis='y',color=GRID,lw=0.6)
for s in ('top','right'): ax.spines[s].set_visible(False)
ax.legend(frameon=False,fontsize=7.6,loc='upper left',labelcolor=INK)
ax.set_title('c   How it is measured (schematic)',loc='left',fontsize=10.5,fontweight='bold',color=INK,x=-0.17)
fig.text(0.03,0.035,'A first-order chain gives a straight line: each extra minute adds the same amount.\nRaw curves rise even without memory (estimator bias), so only the gap to the\nmatched simulation counts.',fontsize=7.8,color=INK2,va='bottom')
# ---------- (d) data: excess per block ----------
O=json.load(open(sqw.OUTD+'p_blockkl.json')); H=json.load(open(sqw.OUTD+'p_blockkl_ddhalves.json'))
REF={(d['g'],d['fly']):d['ref'] for d in O}
def E_of(src,g,p,refkey,nul):
    out=[]
    for d in src:
        if d['g']!=g or REF[(g,d['fly'])][refkey] is None: continue
        s=np.array([x[p]['D3'] for x in d['sims'][nul]]); out.append(d['obs'][p]['D3']-np.median(s))
    return np.array(out)
rows=[('CS LD, light',O,'CS',0,'light','A',False),('CS LD, dark',O,'CS',1,'dark','A',False),
      ('CS DD, ZT 0–12',H,'CSDD',0,'whole','A2',False),('CS DD, ZT 12–24',H,'CSDD',1,'whole','A2',False),
      ('per0 LD, light',O,'per0LD',0,'light','A',False),('per0 LD, dark',O,'per0LD',1,'dark','A',False),
      ('per0 DD',O,'per0DD',0,'whole','A',True),
      ('fumin, light',O,'fumin',0,'light','A',False),('fumin, dark',O,'fumin',1,'dark','A',True),
      ('sss, light',O,'sss',0,'light','A',False),('sss, dark',O,'sss',1,'dark','A',False)]
ax=fig.add_axes([0.66,0.17,0.31,0.30]); rng=np.random.default_rng(1)
for k,(lab,src,g,p,rk,nul,q) in enumerate(rows):
    e=E_of(src,g,p,rk,nul); m=np.median(e); b=[np.median(rng.choice(e,len(e))) for _ in range(2000)]; lo,hi=np.percentile(b,[2.5,97.5])
    y=len(rows)-1-k; ax.plot([lo,hi],[y,y],color=INK2,lw=1.6,solid_capstyle='round')
    ax.plot(m,y,'o',ms=8,mfc=COL['W'] if q else 'white',mec=COL['W'] if q else INK2,mew=1.8)
ax.axvline(0,color=INK2,lw=0.8)
ax.set_yticks(range(len(rows))); ax.set_yticklabels([r[0] for r in rows][::-1],fontsize=8.5)
ax.set_xlabel('extra irreversibility from the 3-min memory, E\n(k$_B$ per hour; median over flies, 95% CI)')
ax.grid(axis='x',color=GRID,lw=0.6)
for s in ('top','right'): ax.spines[s].set_visible(False)
ax.plot([],[],'o',ms=7,mfc=COL['W'],mec=COL['W'],label='passes the pre-registered rule'); ax.plot([],[],'o',ms=7,mfc='white',mec=INK2,mew=1.6,label='does not pass')
ax.legend(frameon=False,fontsize=7.6,loc='lower right',bbox_to_anchor=(1.03,1.0),ncol=2,labelcolor=INK,handletextpad=0.3,columnspacing=1.0)
ax.set_title('d   What we found (20p, data)',loc='left',fontsize=10.5,fontweight='bold',color=INK,x=-0.40,pad=22)
fig.text(0.53,0.035,'fumin dark: the excess sits in windows containing sleep (5-min criterion not excluded).\nper0 DD: present in wake-only windows. CS LD dark: CI excludes 0 but too few single flies (4%) pass.',fontsize=7.8,color=INK2,va='bottom')
os.makedirs('figures',exist_ok=True)
for ext in ('png','pdf'): fig.savefig(f'figures/fig18_memory_dissipation.{ext}',dpi=200,facecolor='white')
print('saved')
