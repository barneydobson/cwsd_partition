# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 17:15:14 2021

@author: bdobson
"""
import os
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib as mpl
from numpy.random import random
import numpy as np

#Compare performance across all partitioning algorithms and number of compartments
#Compare performance across varying timestep size and number of compartments

root = os.path.join("C:\\", "Users", "bdobson", "Documents", "GitHub", "cwsd_sewer","data")
catchment = "cranbrook"
driver = 'GeoJSON'


clusters_root = os.path.join(root,catchment,"results","2021-05-21")

detail_fid = os.path.join(clusters_root, 
                          'repetition_17',
                          'cluster_L_scn_n_48',
                          'sim_dt_60_s',
                          )

plot_folder = os.path.join(clusters_root,"plots")
if not os.path.exists(plot_folder ):
    os.mkdir(plot_folder )

results_arcs = pd.read_csv(os.path.join(clusters_root, "comparison_arcs_nse.csv"))
results_nodes = pd.read_csv(os.path.join(clusters_root, "comparison_nodes_nse.csv"))


def clustercomp2(results,var = 'nse'):
    dt = 60

    ncols = 3
    nrows = 3
    f, axs = plt.subplots(nrows,ncols, figsize=(8,8))
    if 'arc' in results.columns:
        type_ = 'arc'
    else:
        type_ = 'node_id'
    
    results_ = results.copy()
    results_['method_type'] = ['_'.join(x[1:-2]) for x in results_.method.str.split('_')]
    results_ = results_.loc[results_.dt == dt]

    highlight = results_.groupby(type_).median()[var].idxmin()        
    highlight2 = results_.groupby(type_).median()[var].idxmax()        
    print(highlight)
    print(highlight2)
    ind = (results_.method_type == 'louv') & (results_.N < 400) & (results_.N > 300)
    ind = ind | (results_.method_type != 'louv')
    results_ = results_.loc[ind]
    ymin = -25
    yt = [ymin,-1,0,1]
    nb = 5
    
    N = np.logspace(np.log10(results_.N.min()),np.log10(results_.N.max()),nb).astype(int)
    results_['N_plot'] = [N[abs(x-N).argmin()] for x in results_.N]
    
    xt = np.round(np.logspace(np.log10(N[0]*1.05),np.log10(N[-1]*0.95),4)).astype(int)
    
    l=0
    for (idx_clus, grp), ax in zip(results_.groupby('method_type'), axs.reshape(-1)):
        
        for (idx_rain, grp_rain), marker in zip(grp.groupby('rain'), ['x','o']):
            ind1 = grp_rain[type_] == highlight
            ind2 = grp_rain[type_] == highlight2
            ind = ind1 | ind2
            ax.scatter(grp_rain.loc[~ind].N, grp_rain.loc[~ind,var], s=10, marker=marker, color='k')
            ax.scatter(grp_rain.loc[ind1].N, grp_rain.loc[ind1,var], s=10, marker=marker, color='r')
            ax.scatter(grp_rain.loc[ind2].N, grp_rain.loc[ind2,var], s=10, marker=marker, color='c')
        
        
        if idx_clus == 'louv':
            ax.plot([N[0],N[-1]], [grp[var].median(), grp[var].median()],color='b',ls='--')
        else:
            grp_ = grp.groupby('N_plot').median().reset_index()
            ax.plot(grp_.N_plot, grp_[var],color='b',ls='--')
        ax.grid(ls='--')
        ax.set_xscale('log')
        ax.set_yscale('symlog')
        if var == 'nse':
            ax.set_ylim([-20,1])
        else:
            ax.set_ylim([-1,30])
        ax.set_xlim([results.N.min(),results.N.max()])
        ax.set_title(idx_clus)
        ax.set_xticks(xt)
        if int(l / nrows) == (ncols - 1):
            ax.set_xlabel('N compartments')
            ax.set_xticklabels(xt)
        else:
            ax.set_xticklabels([])
        ax.set_yticks(yt)
        if int(l % ncols) == 0:
            ax.set_ylabel(var)
            ax.set_yticklabels(yt)
        else:
            ax.set_yticklabels([])
        l+=1
    f.tight_layout()
    return f
f = clustercomp2(results_arcs, var='nse').savefig(os.path.join(plot_folder, 'cluster_comparison_arcs.svg'))
f = clustercomp2(results_nodes, var= 'nse').savefig(os.path.join(plot_folder, 'cluster_comparison_nodes.svg'))

def dtplot2(results, var = 'nse'):
    f, ax = plt.subplots()
    jitter = 0.1
    cmap = plt.get_cmap('cool')
    
    results_ = results.copy()
    norm = mpl.colors.LogNorm(vmin=results_.N.min(),vmax=results_.N.max())

    nb = 5
    
    N = np.logspace(np.log10(results_.N.min()),np.log10(results_.N.max()),nb).astype(int)
    results_['N_plot'] = [N[abs(x-N).argmin()] for x in results_.N]
    
    for (idx_rain, grp_rain), marker in zip(results_.groupby('rain'), ['o','x']):
        
        for idx_N, grp_N in grp_rain.groupby('N'):
            color = cmap(norm(idx_N))

            ax.plot(grp_N.dt + grp_N.dt * (0.5 - random(grp_N.dt.shape))*jitter, 
                       grp_N[var], 
                       linestyle = 'none', 
                       marker=marker,
                       markersize=3,
                       linewidth=1,
                       markeredgecolor=color,
                       markerfacecolor='none')

        
            
        ax.set_xscale('log')
        if var == 'nse':
            ax.set_ylim([-1,1])
        else:
            ax.set_yscale('symlog')
            ax.set_ylim([-1,30])
    ax.set_ylabel(var)
    ax.set_xlabel('Timestep (seconds)')
    grp_ = results_.groupby(['N_plot','dt']).median().reset_index()
    grp_[var] = grp_[var].clip(lower=-3)
    for N_, group in grp_.groupby('N_plot'):
        color = cmap(norm(N_))
        ax.plot(group.dt, group[var],color=color,ls='--')

    cmaplist = [cmap(norm(x)) for x in N]
    my_cmap = mpl.colors.LinearSegmentedColormap.from_list('Custom cmap', cmaplist, len(N))
    bounds =  np.insert(N+ 1, 0,0)
    my_norm = mpl.colors.BoundaryNorm(bounds, len(bounds))

    cbar = plt.colorbar(mpl.cm.ScalarMappable(norm=my_norm, cmap=my_cmap),
                        ticks = bounds[0:nb] + np.diff(bounds)/2,
                        boundaries = bounds,
                        )
    cbar.set_ticklabels(N)
    
    cbar.ax.set_title('N compartments')
    f.tight_layout()
    return f
dtplot2(results_arcs.loc[results_arcs.method.str.contains('L_L_sc')]).savefig(os.path.join(plot_folder, "timestep_compartment_interaction_arcs.svg"))
dtplot2(results_nodes.loc[results_nodes.method.str.contains('L_L_sc')]).savefig(os.path.join(plot_folder, "timestep_compartment_interaction_nodes.svg"))
