#!/usr/bin/env python3
"""Fixed-scale inspection figures; no data-selected flow or rescaled heatmap colors."""
from pathlib import Path
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

OUT=Path(__file__).resolve().parent
DIAMETERS=[128,235,285,500]
FLOWS=[1,3,5,7,10]
plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False})
read=lambda name:pd.read_csv(OUT/name,float_precision='round_trip')


def save(fig,name):
    fig.savefig(OUT/(name+'.png'),dpi=160,bbox_inches='tight')
    fig.savefig(OUT/(name+'.pdf'),bbox_inches='tight')
    plt.close(fig)


def heat(ax,g,value='rho',coord='absolute',sign=False,title=''):
    matrix=g.pivot(index='source_bin',columns='tail_bin',values=value).sort_index().to_numpy()
    hi=500 if coord=='absolute' else 1
    image=ax.imshow(matrix,aspect='auto',origin='upper',extent=[0,hi,1,0],
        cmap='YlGnBu' if sign else 'RdBu_r',vmin=0 if sign else -1,vmax=5 if sign else 1)
    ax.set_xlabel('Tail depth (um)' if coord=='absolute' else 'Tail depth / diameter')
    ax.set_ylabel('Source normalized depth u');ax.set_title(title)
    ax.set_yticks([0,1/3,2/3,1],['0','1/3','2/3','1'])
    if sign:
        n,m=matrix.shape
        for j in range(n):
            for k in range(m):
                ax.text((k+.5)*hi/m,(j+.5)/n,str(int(matrix[j,k])),ha='center',va='center',fontsize=6,
                        color='white' if matrix[j,k]>=4 else '#222222')
    return image


def main():
    absolute=read('diameter_absolute_consensus_map.csv')
    normalized=read('diameter_normalized_consensus_map.csv')
    maps=read('coupling_maps_all_resolutions_modes.csv.gz')
    for num,data,coord,sign,title in [(1,absolute,'absolute',False,'51-frame mean-mean flow consensus'),
                                   (2,absolute,'absolute',True,'Number positive out of five common flows'),
                                   (4,normalized,'normalized',False,'Direct normalized-tail flow consensus')]:
        fig,axes=plt.subplots(2,2,figsize=(12,7),layout='constrained')
        for ax,d in zip(axes.flat,DIAMETERS):
            im=heat(ax,data[data.diameter_um.eq(d)],'number_positive' if sign else 'rho',coord,sign,f'D{d}')
        fig.colorbar(im,ax=axes.ravel().tolist(),label='Positive volumes / 5' if sign else 'Median descriptive Spearman',shrink=.8)
        fig.suptitle(f'Figure {num}. {title}');save(fig,f'figure{num}_{coord}_'+('positive_consistency' if sign else 'consensus'))
    # Every volume gets raw mean, primary residual mean, and secondary residual Q maps.
    for d in DIAMETERS:
        for scope,flows in [('common_grid',FLOWS)]+([('d500_extension',[2,9,12])] if d==500 else []):
            fig,axes=plt.subplots(3,len(flows),figsize=(3.4*len(flows),8),layout='constrained')
            for col,flow in enumerate(flows):
                for row,(mode,metric,label) in enumerate([('raw','mean','Raw mean'),('detrended_51','mean','51-frame mean'),('detrended_51','q','51-frame Q (secondary)')]):
                    g=maps[maps.diameter_um.eq(d)&maps.flow_mm_s.eq(flow)&maps.coordinate.eq('absolute')&maps.source_bins.eq(6)&maps['mode'].eq(mode)&maps.metric.eq(metric)]
                    im=heat(axes[row,col],g,title=f'{flow} mm/s; {label}')
            fig.colorbar(im,ax=axes.ravel().tolist(),label='Descriptive Spearman; fixed [-1,1]',shrink=.75)
            fig.suptitle(f'Figure 3. D{d}; {scope}');save(fig,f'figure3_D{d}_{scope}_volume_atlas')
    flow=read('flow_map_similarity.csv')
    fig,axes=plt.subplots(1,2,figsize=(11,4.5),layout='constrained')
    for ax,coord in zip(axes,['absolute','normalized']):
        for j,d in enumerate(DIAMETERS):
            g=flow[flow.diameter_um.eq(d)&flow.coordinate.eq(coord)].sort_values(['flow_a','flow_b'])
            ax.scatter(j+np.linspace(-.18,.18,10),g.map_spearman,s=22)
            ax.plot([j-.23,j+.23],[g.map_spearman.median()]*2,color='black',lw=2)
        ax.set_xticks(range(4),[f'D{d}' for d in DIAMETERS]);ax.set_ylim(-1,1);ax.axhline(0,color='gray',lw=.5)
        ax.set_ylabel('Spearman map similarity');ax.set_title(coord+'; ten flow pairs per diameter')
    fig.suptitle('Figure 5. Within-diameter flow-map shape comparison');save(fig,'figure5_flow_map_similarity')
    region=read('volume_high_coupling_region_summary.csv')
    fig,axes=plt.subplots(1,4,figsize=(15,4.5),layout='constrained')
    for ax,d in zip(axes,DIAMETERS):
        g=region[region.diameter_um.eq(d)&region.coordinate.eq('absolute')&region.grid_scope.eq('common_grid')].sort_values('flow_mm_s')
        for color,r in zip(plt.cm.tab10.colors,g.itertuples()):
            ax.errorbar(r.u_high_region_median,r.tail_high_region_median_um,
                xerr=[[r.u_high_region_median-r.u_high_region_lo],[r.u_high_region_hi-r.u_high_region_median]],
                yerr=[[r.tail_high_region_median_um-r.tail_high_region_lo_um],[r.tail_high_region_hi_um-r.tail_high_region_median_um]],
                color=color,marker='o',ms=5,lw=.7,alpha=.75,label=f'{int(r.flow_mm_s)} mm/s')
            ax.scatter(r.u_peak,r.tail_depth_peak_um,marker='x',s=40,color=color)
        ax.set_xlim(0,1);ax.set_ylim(0,500);ax.set_xlabel('Source depth u');ax.set_ylabel('Tail depth (um)');ax.set_title(f'D{d}');ax.legend(fontsize=7)
    fig.suptitle('Figure 6. Circle: high-region median / full bin support; x: maximum cell');save(fig,'figure6_peak_high_region_locations')
    cross=read('absolute_vs_normalized_similarity.csv')
    fig,ax=plt.subplots(figsize=(10,5),layout='constrained')
    for j,r in enumerate(cross.itertuples()):
        ax.plot([j,j],[r.map_spearman_absolute,r.map_spearman_normalized],color='gray',lw=1)
    ax.scatter(range(6),cross.map_spearman_absolute,label='Absolute 6x20',marker='o')
    ax.scatter(range(6),cross.map_spearman_normalized,label='Normalized 6x5',marker='s')
    ax.set_xticks(range(6),[f'D{int(r.diameter_a)} / D{int(r.diameter_b)}' for r in cross.itertuples()],rotation=20)
    ax.set_ylim(-1,1);ax.axhline(0,color='gray',lw=.5);ax.set_ylabel('Consensus-map Spearman similarity');ax.legend()
    ax.set_title('Figure 7. Absolute vs direct normalized-tail coordinates (different support and resolution)')
    save(fig,'figure7_cross_diameter_coordinates')
    # The six cell choices below are fixed in code before reading any correlation values.
    for mode in ['raw','detrended_51']:
        fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained')
        for ax,d in zip(axes.flat,DIAMETERS):
            lag=read(f'lag_absolute_D{d}.csv.gz')
            lag=lag[lag.scan_id.eq(f'D{d}_F05_V01')&lag['mode'].eq(mode)]
            for color,s in zip(['#2678b2','#ba6325','#427c4d'],[2,4,6]):
                for t,ls in [(2,'-'),(10,'--')]:
                    g=lag[lag.source_bin.eq(s)&lag.tail_bin.eq(t)].sort_values('lag')
                    ax.plot(g.lag,g.rho,color=color,ls=ls,lw=1,label=f'S{s}, T{t:02d}')
            ax.axvline(0,color='gray',lw=.6);ax.axhline(0,color='gray',lw=.5);ax.set_ylim(-1,1)
            ax.set_xlabel('Lag in original B-scans');ax.set_ylabel('Descriptive Spearman');ax.set_title(f'D{d}, fixed flow=5');ax.legend(ncol=2,fontsize=7)
        fig.suptitle(f'Figure 8. Fixed cells: S2/S4/S6 x T02(25-50um)/T10(225-250um); {mode}')
        save(fig,'figure8_fixed_flow5_lag_'+mode)
    controlled=read('diameter_geometry_control_consensus.csv')
    fig,axes=plt.subplots(3,4,figsize=(15,8),layout='constrained')
    for col,d in enumerate(DIAMETERS):
        im=heat(axes[0,col],absolute[absolute.diameter_um.eq(d)],title=f'D{d}; uncontrolled')
        for row,control in enumerate(['area_ztop','x1_ztop'],1):
            g=controlled[controlled.diameter_um.eq(d)&controlled.coordinate.eq('absolute')&controlled.control.eq(control)]
            im=heat(axes[row,col],g,title=control)
    fig.colorbar(im,ax=axes.ravel().tolist(),label='Median descriptive (partial) Spearman',shrink=.75)
    fig.suptitle('Figure 9. Geometry controls; signals and controls use 51-frame residuals');save(fig,'figure9_geometry_control_consensus')
    secondary_q_figure()
    print('Fixed-scale figures complete',flush=True)


def secondary_q_figure():
    con=read('diameter_consensus_all_resolutions_modes.csv.gz')
    con=con[con.source_bins.eq(6)&con['mode'].eq('detrended_51')&con.metric.eq('q')]
    fig,axes=plt.subplots(2,4,figsize=(15,6),layout='constrained')
    for row,coord in enumerate(['absolute','normalized']):
        for col,d in enumerate(DIAMETERS):
            im=heat(axes[row,col],con[con.diameter_um.eq(d)&con.coordinate.eq(coord)],coord=coord,title=f'D{d}; {coord}')
    fig.colorbar(im,ax=axes.ravel().tolist(),label='Median descriptive Q-Q Spearman',shrink=.8)
    fig.suptitle('Figure 10. Secondary Q-Q consensus; includes signal and ROI geometry contributions')
    save(fig,'figure10_secondary_q_consensus')


if __name__=='__main__':
    if '--secondary-only' in sys.argv:secondary_q_figure()
    else:main()
