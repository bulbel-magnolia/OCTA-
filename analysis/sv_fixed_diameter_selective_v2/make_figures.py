from common import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
plt.rcParams.update({'font.size':9,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
FIG=OUT/'figures';FIG.mkdir(exist_ok=True)
AUDIT=[]
def save(fig,name,tables,panels):
    fig.savefig(FIG/(name+'.png'),dpi=170,bbox_inches='tight');fig.savefig(FIG/(name+'.pdf'),bbox_inches='tight');plt.close(fig)
    AUDIT.append(dict(figure='figures/'+name+'.png',source_tables=tables,panels=panels,status='generated_from_saved_tables'))
def heat(ax,g,value='rho',excess=False):
    a=g.pivot(index='source_bin',columns='tail_bin',values=value).reindex(index=range(1,7),columns=range(1,21)).to_numpy()
    assert a.shape==(6,20)
    im=ax.imshow(a,origin='upper',extent=[0,500,1,0],aspect='auto',cmap='RdBu_r',vmin=-2 if excess else -1,vmax=2 if excess else 1)
    np.testing.assert_allclose(np.asarray(im.get_array()),a,equal_nan=True);assert im.get_clim()==((-2,2) if excess else (-1,1))
    ax.add_patch(Rectangle((0,1/3),100,2/3,fill=False,edgecolor='black',lw=.8));ax.set(xlabel='Tail band depth (um)',ylabel='Source u',xticks=[0,100,250,500],yticks=[0,1/3,2/3,1]);return im
def main():
    maps=read(OUT/'depth_maps_raw51.csv');controls=read(OUT/'pseudo_maps.csv');ex=read(OUT/'fixed_region_control_cells.csv');cons=read(OUT/'diameter_consensus_maps.csv');lag=read(OUT/'representative_lag.csv.gz')
    for d in [128,235,285]:
        fig,axes=plt.subplots(5,2,figsize=(11,13),layout='constrained')
        for i,f in enumerate([1,3,5,7,10]):
            for j,mode in enumerate(['raw','detrended_51']):
                g=maps[(maps.diameter_um==d)&(maps.flow_mm_s==f)&(maps['mode']==mode)];assert len(g)==120;im=heat(axes[i,j],g);axes[i,j].set_title(f'D{d}, flow {f} | {mode}; +{(g.rho>0).sum()} / -{(g.rho<0).sum()} / NaN {g.rho.isna().sum()}')
        fig.colorbar(im,ax=axes,label='Spearman rho (mean-mean)',shrink=.6);save(fig,f'depth_atlas_D{d}','depth_maps_raw51.csv',10)
        fig,axes=plt.subplots(5,4,figsize=(18,13),layout='constrained')
        for i,f in enumerate([1,3,5,7,10]):
            for j,side in enumerate(['real_matched','left','right']):
                g=controls[(controls.diameter_um==d)&(controls.flow_mm_s==f)&(controls.side==side)];im=heat(axes[i,j],g);axes[i,j].set_title(f'F{f} {side} | 51')
            g=ex[(ex.diameter_um==d)&(ex.flow_mm_s==f)];im2=heat(axes[i,3],g,'excess_cell',True);axes[i,3].set_title(f'F{f} real - (L+R)/2')
        fig.colorbar(im,ax=axes[:,:3],label='rho',shrink=.6);fig.colorbar(im2,ax=axes[:,3],label='Diagnostic excess',shrink=.6);save(fig,f'control_atlas_D{d}','pseudo_maps.csv; fixed_region_control_cells.csv',20)
        fig,axes=plt.subplots(5,3,figsize=(13,12),sharex=True,sharey=True,layout='constrained')
        for i,f in enumerate([1,3,5,7,10]):
            for j,pair in enumerate(REP):
                ax=axes[i,j]
                for mode,color in [('raw','#999999'),('detrended_51','#176b93')]:
                    a=lag[(lag.diameter_um==d)&(lag.flow_mm_s==f)&(lag.pair==pair)&lag['mode'].eq(mode)].sort_values('lag');assert len(a)==101
                    ax.plot(a.lag,a.rho,color=color,label=mode)
                ax.axvline(0,c='black',lw=.6);ax.axhline(0,c='gray',lw=.5);ax.set(title=f'F{f} {pair}',ylim=(-1,1),xlabel='Lag (original frames)',ylabel='Spearman rho')
        axes[0,0].legend(fontsize=7);save(fig,f'lag_atlas_D{d}','representative_lag.csv.gz',15)
    fig,axes=plt.subplots(3,2,figsize=(11,8),layout='constrained')
    for i,d in enumerate([128,235,285]):
        for j,mode in enumerate(['raw','detrended_51']):
            g=cons[(cons.diameter_um==d)&cons['mode'].eq(mode)];im=heat(axes[i,j],g,'median');axes[i,j].set_title(f'D{d} {mode}: median of 5 flow coefficients')
    fig.colorbar(im,ax=axes,label='Median rho; no pooled frames',shrink=.6);save(fig,'consensus','diameter_consensus_maps.csv',6)
    volume=read(OUT/'volume_metrics.csv');fig,axes=plt.subplots(1,5,figsize=(16,3.6),layout='constrained')
    for ax,metric in zip(axes,['source','tail100','tail500','RI100','RI500']):
        for d,color in [(128,'#176b93'),(235,'#009e73'),(285,'#d55e00')]:
            a=volume[(volume.diameter_um==d)&volume.metric.eq(metric)].sort_values('flow_mm_s');ax.plot(a.flow_mm_s,a['median'],'o-',color=color,label=f'D{d}')
        ax.set(title=metric,xlabel='Flow (mm/s)',ylabel='Tail/source' if metric.startswith('RI') else 'Raw SV instrument units');ax.set_ylim(bottom=0)
    axes[0].legend();save(fig,'matched_flow_metrics','volume_metrics.csv',5)
    depth=read(OUT/'raw_depth_profiles.csv');local=read(OUT/'local_contrast_profiles.csv');fig,axes=plt.subplots(3,3,figsize=(13,10),layout='constrained')
    for i,d in enumerate([128,235,285]):
        for j,metric in enumerate(['tail_mean','RI','z_local']):
            ax=axes[i,j]
            for f in [1,3,5,7,10]:
                if metric=='z_local':
                    a=local[(local.diameter_um==d)&local.flow_mm_s.eq(f)&local.metric.eq(metric)&local.region.str.fullmatch('t[0-9]{2}')].copy();a['tail_mid_um']=(a.region.str[1:].astype(int)-.5)*25
                else:a=depth[(depth.diameter_um==d)&depth.flow_mm_s.eq(f)&depth.metric.eq(metric)]
                a=a.sort_values('tail_mid_um');assert len(a)==20;ax.plot(a.tail_mid_um,a['median'],'o-',ms=2,lw=1,label=f'F{f}')
            ax.axhline(0,c='gray',lw=.6);ax.set(title=f'D{d}: {metric}',xlabel='Tail band midpoint (um)',ylabel='Volume median '+metric,xlim=(0,500))
        axes[i,0].legend(ncol=3,fontsize=7)
    save(fig,'depth_profiles','raw_depth_profiles.csv; local_contrast_profiles.csv',9)
    e=read(OUT/'fixed_region_control_excess.csv');fig,axes=plt.subplots(1,3,figsize=(13,4),layout='constrained')
    for ax,d in zip(axes,[128,235,285]):
        g=e[e.diameter_um==d].sort_values('flow_mm_s')
        for field,label in [('real_matched_rho_value','Real C16'),('left_rho_value','L C16'),('right_rho_value','R C16'),('excess_cell_value','C16 cellwise excess')]:ax.plot(g.flow_mm_s,g[field],'o-',label=label)
        ax.axhline(0,color='gray',lw=.7);ax.set(title=f'D{d}',xlabel='Flow (mm/s)',ylabel='Signed coefficient summary',ylim=(-1,1))
    axes[0].legend(fontsize=8);save(fig,'fixed_region_controls','fixed_region_control_excess.csv',3)
    # Position and detrend changes are paired coefficient differences, not CIs.
    pp=read(OUT/'position_sensitivity_pairs.csv');fig,axes=plt.subplots(3,3,figsize=(13,10),layout='constrained')
    for i,d in enumerate([128,235,285]):
        for j,rep in enumerate(REP):
            ax=axes[i,j]
            for name in ['baseline','xp','xm','zp','zm']:
                g=pp[(pp.diameter_um==d)&pp.pair.eq(rep)&pp.offset.eq(name)].sort_values('flow_mm_s');ax.plot(g.flow_mm_s,g.rho0,'o-',label=name,lw=1)
            ax.set(title=f'D{d} {rep}',ylim=(-1,1),xlabel='Flow (mm/s)',ylabel='51-frame rho, common support')
        axes[i,0].legend(ncol=3,fontsize=7)
    save(fig,'position_sensitivity','position_sensitivity_pairs.csv',9)
    csv('figure_manifest.csv',AUDIT);js('figure_validation.json',dict(status='generated_and_numeric_mapping_passed',figures=len(AUDIT),all_main_map_cells_displayed=True,all_map_limits_rho=[-1,1],all_excess_limits=[-2,2],physical_depth_last_band=[475,500],visual_review='pending'))
if __name__=='__main__':main()
