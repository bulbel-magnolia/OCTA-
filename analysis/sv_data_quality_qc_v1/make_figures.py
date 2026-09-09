"""Fixed cohort figures: independent volumes are the displayed replicate units."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from derive_qc import OUT,read
DS=[128,235,285,500];COLORS=['#0072B2','#009E73','#D55E00','#CC79A7']
plt.rcParams.update({'font.size':9,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
def save(fig,n):
    fig.savefig(OUT/(n+'.png'),dpi=180,bbox_inches='tight');fig.savefig(OUT/(n+'.pdf'),bbox_inches='tight');plt.close(fig)
def zero_limits(a):
    a=np.asarray(a,float);a=a[np.isfinite(a)];lo=min(0,a.min());hi=max(0,a.max());pad=(hi-lo)*.06 or 1
    return lo-pad if lo<0 else 0,hi+pad
def main():
    v=read(OUT/'volume_qc_summary.csv');v=v[v.grid_scope.eq('common_grid')]
    r=read(OUT/'volume_region_detectability_summary.csv');r=r[r.grid_scope.eq('common_grid')]
    n=read(OUT/'negative_control_volume_summary.csv');n=n[n.grid_scope.eq('common_grid')]
    fig,ax=plt.subplots(figsize=(7,4))
    for d,c in zip(DS,COLORS):
        a=v[v.diameter_um.eq(d)];ax.plot(a.flow_mm_s,a.source_z,'o-',color=c,label=f'D{d}')
    ax.set(xlabel='Flow (mm/s)',ylabel='Whole-source robust Z',ylim=zero_limits(v.source_z));ax.legend(ncol=4);save(fig,'figure1_source_detectability')
    t=r[r.region.str.fullmatch('t[0-9]{2}')].copy();t['depth']=(t.region.str[1:].astype(int)-.5)*25
    fig,axes=plt.subplots(2,2,figsize=(10,7),sharex=True,sharey=True)
    for d,c,ax in zip(DS,COLORS,axes.flat):
        a=t[t.diameter_um.eq(d)]
        for f,b in a.groupby('flow_mm_s'):ax.plot(b.depth,b.z_median,color=c,alpha=.35,lw=.9,label=f'F{f}')
        m=a.groupby('depth').z_median.median();ax.plot(m.index,m.values,color=c,lw=2.5,label='5-volume median');ax.set_title(f'D{d}')
        ax.set(xlabel='Tail depth (um)',ylabel='Volume median robust Z',ylim=zero_limits(t.z_median));ax.axhline(0,color='gray',lw=.6)
    axes[0,0].legend(ncol=3,fontsize=7);save(fig,'figure2_tail_depth_detectability')
    fig,axes=plt.subplots(1,5,figsize=(13,3.3),sharey=True)
    regions=['t01','t02','t03','t04','tail100'];labels=['0-25','25-50','50-75','75-100','Pooled 0-100']
    lim=zero_limits(r[r.region.isin(regions)].z_median)
    for name,label,ax in zip(regions,labels,axes):
        for i,(d,c) in enumerate(zip(DS,COLORS)):
            a=r[r.region.eq(name)&r.diameter_um.eq(d)].sort_values('flow_mm_s');ax.scatter(i+np.linspace(-.13,.13,len(a)),a.z_median,color=c,s=20);ax.hlines(a.z_median.median(),i-.25,i+.25,color=c,lw=2)
        ax.set(title=label+' um',xticks=range(4),xticklabels=DS,ylim=lim)
    axes[0].set_ylabel('Volume median robust Z');save(fig,'figure3_proximal_detectability')
    s=r[r.region.str.fullmatch('s[1-6]')].pivot(index='scan_id',columns='region',values='z_median')
    fig,ax=plt.subplots(figsize=(8,7));im=ax.imshow(s.values,aspect='auto',cmap='viridis',vmin=0,vmax=max(1,s.values.max()));ax.set(xticks=range(6),xticklabels=s.columns,yticks=range(len(s)),yticklabels=s.index,xlabel='Normalized source axial bin');fig.colorbar(im,ax=ax,label='Volume median robust Z');save(fig,'figure4_source_depth_detectability')
    fig,axes=plt.subplots(1,3,figsize=(12,3.8))
    for ax,col,label in zip(axes,['background_rcv','background_drift','background_trend_amplitude'],['Background RCV','First-to-last relative drift','51-frame trend amplitude / median']):
        for d,c in zip(DS,COLORS):
            a=v[v.diameter_um.eq(d)];ax.plot(a.flow_mm_s,a[col],'o-',color=c,label=f'D{d}')
        ax.set(xlabel='Flow (mm/s)',ylabel=label,ylim=zero_limits(v[col]))
    axes[0].legend();save(fig,'figure5_background_stability')
    fig,axes=plt.subplots(1,3,figsize=(12,3.8))
    for ax,col in zip(axes,['X4_jitter','X1_jitter','z_top_jitter']):
        for d,c in zip(DS,COLORS):
            a=v[v.diameter_um.eq(d)];ax.plot(a.flow_mm_s,a[col],'o-',color=c,label=f'D{d}')
        ax.set(xlabel='Flow (mm/s)',ylabel='Median absolute consecutive change (pixel)',title=col,ylim=zero_limits(v[col]))
    axes[0].legend();save(fig,'figure6_geometry_jitter')
    fig,axes=plt.subplots(2,2,figsize=(10,7),sharex=True,sharey=True)
    for d,ax in zip(DS,axes.flat):
        a=n[n.diameter_um.eq(d)]
        for col,label,c in [('real_high_median','Real high region','#222222'),('left_high_at_real_region','Left at real region','#0072B2'),('right_high_at_real_region','Right at real region','#D55E00')]:ax.plot(a.flow_mm_s,a[col],'o-',label=label,color=c)
        ax.set(title=f'D{d}',xlabel='Flow (mm/s)',ylabel='51-frame median rho at real high cells',ylim=(-1,1));ax.axhline(0,color='gray',lw=.7)
    axes[0,0].legend(fontsize=8);save(fig,'figure7_real_vs_pseudo')
    e=read(OUT/'control_excess_coupling_summary.csv');e=e[e.grid_scope.eq('common_grid')]
    fig,axes=plt.subplots(4,5,figsize=(16,10),sharex=True,sharey=True,layout='constrained')
    for row,d in enumerate(DS):
        for col,f in enumerate([1,3,5,7,10]):
            a=e[e.diameter_um.eq(d)&e.flow_mm_s.eq(f)].pivot(index='source_bin',columns='tail_bin',values='control_excess_rho')
            im=axes[row,col].imshow(a,aspect='auto',extent=[0,500,1,0],vmin=-2,vmax=2,cmap='RdBu_r');axes[row,col].set_title(f'D{d}, F{f}')
            if col==0:axes[row,col].set_ylabel('Source u')
            if row==3:axes[row,col].set_xlabel('Tail depth (um)')
    fig.colorbar(im,ax=axes,label='Real rho minus median pseudo rho; diagnostic only',shrink=.7);save(fig,'figure8_control_excess_maps')
    predictors=['source_z','proximal_z','background_rcv','z_top_jitter','X1_jitter','valid_fraction','negative_control_rho']
    fig,axes=plt.subplots(4,7,figsize=(18,10),sharey=True)
    for i,(d,c) in enumerate(zip(DS,COLORS)):
        a=v[v.diameter_um.eq(d)]
        for j,p in enumerate(predictors):
            ax=axes[i,j];ax.scatter(a[p],a.real_high_rho,color=c)
            for _,x in a.iterrows():ax.annotate(str(int(x.flow_mm_s)),(x[p],x.real_high_rho),fontsize=6)
            ax.set(xlim=zero_limits(v[p]),ylim=zero_limits(v.real_high_rho))
            if i==3:ax.set_xlabel(p,fontsize=8)
            if j==0:ax.set_ylabel(f'D{d}: real high rho')
    save(fig,'figure9_qc_vs_coupling')
    fig,axes=plt.subplots(2,4,figsize=(13,7))
    for ax,p in zip(axes.flat,predictors+['real_high_rho']):
        for i,(d,c) in enumerate(zip(DS,COLORS)):
            a=v[v.diameter_um.eq(d)].sort_values('flow_mm_s');ax.scatter(i+np.linspace(-.13,.13,len(a)),a[p],color=c);ax.hlines(a[p].median(),i-.24,i+.24,color=c,lw=2)
        ax.set(title=p,xticks=range(4),xticklabels=DS,ylim=zero_limits(v[p]))
    fig.suptitle('Common-flow volumes only; points = scans, bars = five-volume medians');fig.tight_layout();save(fig,'figure10_d500_qc_comparison')

if __name__=='__main__':main()
