# Figures for PRESPEC 20h: dwell times. Fig 9 observed S/wake/Q/W survival with model; 9s hazards with model band;
# Fig 10 hidden H1/H2/H3 decoded survival vs geometric + mean dwell; 10s hidden histograms.
import json,numpy as np,matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import sqw_ref as R
plt.rcParams.update({'font.family':'STIXGeneral','mathtext.fontset':'stix','font.size':8,'axes.linewidth':0.6,'xtick.direction':'in','ytick.direction':'in'})
S=json.load(open('sqw_out/h_summary.json')); INK='#1a1a1a'
COL={'CS':'#2a78d6','fumin':'#eb6834','per0LD':'#1baf7a','per0DD':'#eda100','sss':'#e87ba4'}
NAME={'CS':'CS','fumin':'$fmn$','per0LD':'$per^0$ LD','per0DD':'$per^0$ DD','sss':'$sss$'}
ROWS=[('light',[('CS','light'),('per0LD','light'),('fumin','light'),('sss','light')]),
      ('dark / DD',[('CS','dark'),('per0LD','dark'),('per0DD','whole'),('fumin','dark'),('sss','dark')])]
E=np.array(R.DWELL_EDGES); mid=np.sqrt(E[:-1]*np.maximum(E[1:]-1,E[:-1]))
# ---- Fig 9 ----
fig,axs=plt.subplots(2,4,figsize=(7.0,3.9),sharex=True,sharey=True)
for r,(rn,blocks) in enumerate(ROWS):
    for c,nm in enumerate(['S','wake','Q','W']):
        ax=axs[r,c]
        for g,bn in blocks:
            s=S[f'obs|{g}|{bn}|{nm}']; t=np.array(s['t']); y=np.array(s['S'])
            ax.step(t,np.maximum(y,1e-4),where='post',color=COL[g],lw=1.0)
            ax.plot(np.arange(1,721),np.maximum(s['sim_S_med'],1e-4),color=COL[g],lw=0.7,ls=(0,(2,1.5)))
        ax.set_xscale('log'); ax.set_yscale('log'); ax.set_ylim(1e-3,1.05); ax.set_xlim(1,720)
        if r==0: ax.set_title({'S':'sleep (S)','wake':'wake (Q or W)','Q':'quiet wake (Q)','W':'active wake (W)'}[nm],fontsize=8)
        if c==0: ax.set_ylabel(f'{rn}\nP(run > t)')
        if r==1: ax.set_xlabel('t (min)')
for g in COL: axs[0,3].plot([],[],color=COL[g],lw=1.2,label=NAME[g])
axs[0,3].plot([],[],color=INK,lw=0.7,ls=(0,(2,1.5)),label='fitted model'); axs[0,3].legend(frameon=False,fontsize=6,loc='lower left')
fig.tight_layout(pad=0.4); fig.savefig('figures/fig9_dwell_observed.pdf'); fig.savefig('figures/fig9_dwell_observed.png',dpi=300); plt.close(fig)
# ---- Fig 9s: hazards with model band ----
B=[(g,bn) for _,bl in ROWS for g,bn in bl]
fig,axs=plt.subplots(len(B),4,figsize=(7.0,11),sharex=True)
for r,(g,bn) in enumerate(B):
    for c,nm in enumerate(['S','wake','Q','W']):
        ax=axs[r,c]; s=S[f'obs|{g}|{bn}|{nm}']; h=np.array(s['h']); show=np.array(s['n_risk'])>=10
        ax.fill_between(mid,np.array(s['h_lo']),np.array(s['h_hi']),color='#d9d9d9',lw=0,step=None)
        ax.plot(mid[show],h[show],'o-',color=COL[g],ms=2.5,lw=0.8)
        ax.set_xscale('log'); ax.set_ylim(0,max(0.05,np.nanmax(np.r_[h[show],np.array(s['h_hi'])[show]])*1.1))
        if r==0: ax.set_title({'S':'sleep','wake':'wake','Q':'quiet','W':'active'}[nm],fontsize=8)
        if c==0: ax.set_ylabel(f"{NAME[g]} {bn}\nhazard (min$^{{-1}}$)",fontsize=6.5)
        ax.tick_params(labelsize=6)
for ax in axs[-1]: ax.set_xlabel('t (min)')
fig.tight_layout(pad=0.3); fig.savefig('figures/fig9s_hazard_observed.pdf'); fig.savefig('figures/fig9s_hazard_observed.png',dpi=250); plt.close(fig)
# ---- Fig 10: hidden ----
fig=plt.figure(figsize=(7.0,4.4)); gs=fig.add_gridspec(2,4,width_ratios=[1,1,1,1.15])
for r,(rn,blocks) in enumerate(ROWS):
    for c,nm in enumerate(['H1','H2','H3']):
        ax=fig.add_subplot(gs[r,c])
        for g,bn in blocks:
            s=S[f'hid|{g}|{bn}|{nm}']; ax.step(s['t'],np.maximum(s['S'],1e-4),where='post',color=COL[g],lw=1.0)
            ax.plot(np.arange(0,721),np.maximum(s['geo_S'],1e-4),color=COL[g],lw=0.7,ls=(0,(2,1.5)))
        ax.set_xscale('log'); ax.set_yscale('log'); ax.set_ylim(1e-3,1.05); ax.set_xlim(1,720)
        if r==0: ax.set_title({'H1':'H1 (quiet / regulator)','H2':'H2 (wake)','H3':'H3 (sleep)'}[nm],fontsize=8)
        if c==0: ax.set_ylabel(f'{rn}\nP(dwell > t)')
        if r==1: ax.set_xlabel('t (min)')
ax=fig.add_subplot(gs[:,3]); xs=np.arange(3); k=0
allb=[(g,bn) for _,bl in ROWS for g,bn in bl]
for j,(g,bn) in enumerate(allb):
    vals=[np.median(S[f'hid|{g}|{bn}|{nm}']['mean_dwell_fly']) for nm in ['H1','H2','H3']]
    off=(j-len(allb)/2)*0.08
    ax.plot(xs+off,vals,'o',ms=4,mfc=COL[g] if bn!='light' else 'white',mec=COL[g],mew=1)
ax.set_yscale('log'); ax.set_xticks(xs); ax.set_xticklabels(['H1','H2','H3']); ax.set_ylabel('mean dwell 1/(1−$P_{ii}$) (min), fly median')
ax.plot([],[],'o',mfc='white',mec=INK,ms=4,label='light'); ax.plot([],[],'o',mfc=INK,mec=INK,ms=4,label='dark / DD')
for g in COL: ax.plot([],[],'s',color=COL[g],ms=4,label=NAME[g])
ax.legend(frameon=False,fontsize=6,loc='upper left')
fig.text(0.01,0.005,'Solid: decoded posterior paths. Dashed: geometric dwell implied by each fly\'s $P_{ii}$, averaged over flies weighted by their decoded run counts.',fontsize=6,color='#333333')
fig.tight_layout(rect=(0,0.03,1,1),pad=0.4); fig.savefig('figures/fig10_dwell_hidden.pdf'); fig.savefig('figures/fig10_dwell_hidden.png',dpi=300); plt.close(fig)
# ---- Fig 10s: hidden histograms (density per minute) ----
fig,axs=plt.subplots(2,3,figsize=(7.0,4.0),sharex=True)
for r,(rn,blocks) in enumerate(ROWS):
    for c,nm in enumerate(['H1','H2','H3']):
        ax=axs[r,c]
        for g,bn in blocks:
            hst=np.array(S[f'hid|{g}|{bn}|{nm}']['hist']); dens=hst/hst.sum()/np.diff(E)
            ax.plot(mid,np.maximum(dens,1e-6),'o-',ms=2.5,lw=0.8,color=COL[g])
        ax.set_xscale('log'); ax.set_yscale('log'); ax.set_ylim(1e-5,1)
        if r==0: ax.set_title(nm,fontsize=8)
        if c==0: ax.set_ylabel(f'{rn}\ndensity (min$^{{-1}}$)')
        if r==1: ax.set_xlabel('dwell (min)')
fig.tight_layout(pad=0.4); fig.savefig('figures/fig10s_dwell_hidden_hist.pdf'); fig.savefig('figures/fig10s_dwell_hidden_hist.png',dpi=250); plt.close(fig)
print('ok')
