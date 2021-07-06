# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25 16:14:45 2021

@author: bdobson
"""

import os
import pandas as pd
import geopandas as gpd
from matplotlib import pyplot as plt


root = os.path.join("C:\\", "Users", "bdobson", "Documents", "GitHub", "cwsd_sewer","data")
catchment = "cranbrook"
cluster = 'cluster_Louv_266'
rain = "august"
dt = "sim_dt_30_s"
cluster_root = os.path.join(root,catchment,"results","2021-03-02",cluster,dt)
results_root = os.path.join(cluster_root, rain)
plots_root = os.path.join(root, catchment, "results", "2021-03-02", "plots")
info_fid = os.path.join(results_root, "highfid_flows.gzip")
flow_fid = os.path.join(results_root, "flows.gzip")
dep_fid = os.path.join(results_root, "depths.gzip")

info_df = pd.read_parquet(info_fid).set_index('time')
info_df.index = pd.to_datetime(info_df.index)

flow_df = pd.read_parquet(flow_fid)
flow_df.index = pd.to_datetime(flow_df.index)

edges_gdf = gpd.read_file(os.path.join(cluster_root, "compartment_edges.geojson"))
nodes_gdf = gpd.read_file(os.path.join(cluster_root, "compartment_nodes.geojson"))

info_fid = os.path.join(results_root, "highfid_nodes.gzip")
node_fid = os.path.join(results_root, "storages.gzip")

infon_df = pd.read_parquet(info_fid).set_index('time')
infon_df.index = pd.to_datetime(infon_df.index)
# infon_df['volume'] = infon_df['volume'] + infon_df['floodvol']
node_df = pd.read_parquet(node_fid)
node_df.index = pd.to_datetime(node_df.index)
node_df.node = node_df.node.str.replace('-sewer','').astype(int)

dep_df = pd.read_parquet(dep_fid)
dep_df.index = pd.to_datetime(dep_df.index)
dep_df.node = dep_df.node.str.replace('-sewer','').astype(int)

flood_fid = os.path.join(results_root, "flood_vol.gzip")
flood_df = pd.read_parquet(flood_fid)
flood_df.index = pd.to_datetime(flood_df.index)
flood_df.node = flood_df.node.str.replace('-land','').astype(int)
def plot_arc(names, dr):
    f, axs = plt.subplots(len(names),1,figsize=(5,7.5))
    # for i, name in enumerate(names):
    #     ax = axs[i,1]
    for name, ax in zip(names, axs):
        nse = pd.merge(flow_df.loc[flow_df.arc == name, 'flow_out'].rolling('300s').mean().reset_index(), info_df.loc[info_df.info_id == name,'val'].reset_index(), left_on = 'date', right_on = 'time')
        nse = nse.dropna()
        nse = 1 - sum((nse.flow_out - nse.val)**2) / sum((nse.val - nse.val.mean())**2)
        print(nse)
        flow_df.loc[(flow_df.arc == name) & (flow_df.index.isin(pd.date_range(dr[0],dr[1],freq='60s'))), ['flow_out']].plot(color = ['c'], ax=ax)
        # flow_df.loc[(flow_df.arc == name) & (flow_df.index.isin(pd.date_range(dr[0],dr[1],freq='60s'))), ['flow_out']].plot(color = ['c'], ax=ax)
        
        info_ss = info_df.loc[(info_df.info_id == name) & (info_df.index.isin(pd.date_range(dr[0],dr[1],freq='60s'))), 'val']
        info_ss = info_ss.resample('60s').mean().fillna(0)
        info_ss.plot(color = 'k', ax=ax,marker='.',linestyle='',markersize=1)
        ax.set_xlim([dr[0],dr[1]])
        
        if name != names[-1]:
                    ax.set_xlabel('')
                    ax.set_xticklabels([])
                    ax.set_xticks([])
        else:
            ax.set_xlabel('Time (min)')
        
        # ax = axs[i,0]
        # flow_df.loc[flow_df.arc == name, ['flow_out']].rolling('300s').mean().plot(color = ['c'], ax=ax)
        
        # info_df.loc[info_df.info_id == name, 'val'].resample('60s').mean().fillna(0).plot(color = 'k', ax=ax,marker='.',linestyle='',markersize=1)
        
        
        # mm = flow_df.loc[(flow_df.arc == name) & (flow_df.index.isin(pd.date_range(dr[0],dr[1],freq='60s'))), 'flow_out'].max()
        ax.set_ylabel('Flow (m3/min)')
        if name != names[-1]:
            ax.set_xlabel('')
            ax.set_xticklabels([])
            ax.set_xticks([])
        else:
            ax.set_xlabel('Time (min)')
        ax.legend(['CWSD','InfoWorks'])
        
        # ax.set_ylim([0, mm*1.75])
        
    f.tight_layout()
    
    return f


def plot_node(names, df, info_lab,dr):
    f, axs = plt.subplots(len(names),1,figsize=(5,7.5))
    for name, ax in zip(names,axs):
        nse = pd.merge(df.loc[df.node == name, 'val'].reset_index(), infon_df.loc[infon_df[cluster] == name,'volume'].reset_index(), left_on = 'date', right_on = 'time')
        nse = nse.dropna()
        nse = 1 - sum((nse.val - nse.volume)**2) / sum((nse.volume- nse.volume.mean())**2)
        print(nse)
        
        df.loc[(df.node == name) & (df.index.isin(pd.date_range(dr[0],dr[1],freq='60s'))), 'val'].plot(color = 'c', ax=ax)
        
        infon_ss = infon_df.loc[(infon_df[cluster] == name) & (infon_df.index.isin(pd.date_range(dr[0],dr[1],freq='60s'))), info_lab]
        infon_ss = infon_ss.resample('60s').mean().fillna(0)
        infon_ss.plot(color = 'k', ax=ax,marker='.',linestyle='',markersize=1)
        # mm = df.loc[(df.node == name) & (df.index.isin(pd.date_range(dr[0],dr[1],freq='30s'))), 'val'].max()
        ax.set_ylabel('Volume (m3)')
        if name != names[-1]:
                    ax.set_xlabel('')
                    ax.set_xticklabels([])
                    ax.set_xticks([])
        else:
            ax.set_xlabel('Time (min)')

        ax.legend(['CWSD','InfoWorks'])
        ax.set_xlim([dr[0],dr[1]])
        # ax.set_ylim([0, mm*1.3])
    f.tight_layout()
    return f

# f = plot_arc(['node_106.1','node_817.1','node_1439.1'], )
f = plot_arc(['node_106.1','node_817.2','node_1753.1'], [pd.to_datetime('2017-08-07 00:00'), pd.to_datetime('2017-08-14 00:00')])#.savefig(os.path.join(plots_root, "flow_example.svg"))
mdf = pd.merge(node_df.reset_index(), flood_df.reset_index(),on=['date','node']).set_index('date')
mdf['val'] = mdf.val_x + mdf.val_y
f = plot_node([1,20,17], mdf, 'volume', [pd.to_datetime('2017-08-07 00:00'), pd.to_datetime('2017-08-14 00:00')])#.savefig(os.path.join(plots_root, "node_example.svg"))


f = plot_arc(['node_1863.1','node_2477.1','lake_3.3','lake_3.2','lake_3.1'], [pd.to_datetime('2017-08-07 00:00'), pd.to_datetime('2017-08-14 00:00')])
f = plot_node([17,250], mdf, 'volume', [pd.to_datetime('2017-08-07 00:00'), pd.to_datetime('2017-08-14 00:00')])

sum('asd')
f.savefig(os.path.join(results_root, "example_flows.svg"))

mdf = pd.merge(node_df.reset_index(), flood_df.reset_index(),on=['date','node']).set_index('date')
mdf['val'] = mdf.val_x + mdf.val_y
f = plot_node([1,20,17], mdf, 'volume', [pd.to_datetime('2017-08-09 00:00'), pd.to_datetime('2017-08-10 00:00')]).savefig(os.path.join(plots_root, "node_example.svg"))


f = plot_node([15,32], mdf, 'volume', [pd.to_datetime('2017-08-09 00:00'), pd.to_datetime('2017-08-10 00:00')])
f.savefig(os.path.join(results_root, "example_nodes.svg"))


f = plot_node([15,32], dep_df, 'depth', [pd.to_datetime('2017-08-09 00:00'), pd.to_datetime('2017-08-10 00:00')])

#Flooded comparison
fv = infon_df.dropna().groupby(cluster).floodvol.max()

mdf = pd.merge(mdf, nodes_gdf[['storage', 'compart_id']], left_on = 'node', right_on = 'compart_id', how = 'left')
mdf['percent_full'] = mdf['val'] / mdf['storage']

max_val = mdf.groupby('node').val.max()

xy = pd.merge(fv.reset_index(), max_val.reset_index(), left_on = cluster, right_on = 'node')

plt.scatter(xy.val, xy.floodvol)
maxo = xy.max().max()
plt.plot([0,maxo], [0,maxo], ls='--',color='r')
plt.yscale('symlog')
plt.xscale('symlog')
plt.xlabel('Max CWSD compartment volume (m3)')
plt.ylabel('Max ICM flooded volume (m3)')
plt.axis('square')

plt.savefig(os.path.join(plots_root, "flood_comparison.svg"))

