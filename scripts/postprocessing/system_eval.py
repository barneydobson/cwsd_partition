
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  4 15:39:23 2021

@author: bdobson
"""

import os
import pandas as pd
import geopandas as gpd
from matplotlib import pyplot as plt

#Evaluate performance metrics for all nodes/arcs and merge them into geojson 
#for saving

#Plot relationship between node/arc size and NSE 


root = os.path.join("C:\\", "Users", "bdobson", "Documents", "GitHub", "cwsd_sewer","data")
catchment = "cranbrook"
cluster = 'cluster_Louv_266'

for rain, dt in zip(['august','january','august','january'], ['30','60','60','30']):
    driver = 'GeoJSON'
    cluster_root = os.path.join(root,catchment,"results","2021-03-02",cluster,"sim_dt_" + dt + "_s")
    results_root = os.path.join(cluster_root,rain)
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
    
    node_df = pd.read_parquet(node_fid)
    node_df.index = pd.to_datetime(node_df.index)
    node_df.node = node_df.node.str.replace('-sewer','').astype(int)
    
    flood_fid = os.path.join(results_root, "flood_vol.gzip")
    flood_df = pd.read_parquet(flood_fid)
    flood_df.index = pd.to_datetime(flood_df.index)
    flood_df.node = flood_df.node.str.replace('-land','').astype(int)
    
    flow_df['flow'] = flow_df['flow_out']
    
    mdf = pd.merge(node_df.reset_index(), flood_df.reset_index(),on=['date','node']).set_index('date')
    mdf['val'] = mdf.val_x + mdf.val_y
    
    flows = pd.merge(flow_df[['arc','flow']].reset_index(), info_df[['arc','val']].reset_index(), left_on = ['date','arc'], right_on = ['time','arc'], how = 'right').dropna()
    def f(x):
        if x.val.var() > 0:
            obj = 1 - sum((x.flow - x.val)**2)/sum((x.val - x.val.mean())**2)
        else:
            obj = None
        return obj
    nse_flows = flows.groupby('arc').apply(f)
    

    def f(x):
        if x.volume.var() > 0:
            obj = 1 - sum((x.volume - x.val)**2)/sum((x.volume - x.volume.mean())**2)
        else:
            obj = None
        return obj
    nodes = pd.merge(mdf[['node','val']].reset_index(), infon_df[[cluster,'volume']].reset_index(), left_on = ['date','node'], right_on = ['time', cluster]).dropna()
    nodes.loc[nodes.val == 0, 'val'] = 0.001
    nodes.loc[nodes.volume == 0, 'val'] = 0.001
    
    nse_nodes = nodes.groupby('node').apply(f)
    
    nse_comp = pd.merge(flow_df.groupby('arc').max().flow.reset_index(), nse_flows.rename('nse').reset_index(), on='arc')
    f, axs = plt.subplots(1,2,figsize=(8,3.5))
    axs[0].axvspan(xmin=0,xmax=1,color='r',alpha=0.1)
    axs[0].axvspan(xmin=1,xmax=10,color='b',alpha=0.1)
    axs[0].axvspan(xmin=10,xmax=1442,color='g',alpha=0.1)
    axs[0].scatter(nse_comp.flow, nse_comp.nse, s = 10, c='k')
    axs[0].set_xscale('symlog')
    axs[0].set_yscale('symlog')
    axs[0].set_ylim([-10,1])
    axs[0].set_xlabel('Max Pipe Flow (m3/min)')
    axs[0].set_ylabel('NSE')
    
    nse_comp = pd.merge(nodes.groupby('node').max().volume.reset_index(), nse_nodes.rename('nse').reset_index(), on='node')
    
    axs[1].axvspan(xmin=0,xmax=100,color='r',alpha=0.1)
    axs[1].axvspan(xmin=100,xmax=1000,color='b',alpha=0.1)
    axs[1].axvspan(xmin=1000,xmax=50000,color='g',alpha=0.1)
    
    axs[1].scatter(nse_comp.volume, nse_comp.nse, s = 10, c='k')
    axs[1].set_xscale('symlog')
    axs[1].set_yscale('symlog')
    axs[1].set_ylim([-10,1])
    axs[1].set_xlim([3,50000])
    axs[1].set_xlabel('Max Compartment Volume (m3)')
    axs[1].set_ylabel('NSE')
    axs[1].grid(False)
    f.tight_layout()
    f.savefig(os.path.join(results_root,'size_comparison.svg'))
    
    edges_gdf_results = pd.merge(edges_gdf, nse_flows.reset_index(), left_on = 'info_id', right_on = 'arc').rename(columns={0 : 'nse'})
    edges_gdf_results = pd.merge(edges_gdf_results, flow_df.groupby('arc').flow_out.max().reset_index(), how='left')
    edges_gdf_results.to_file(os.path.join(results_root, "flows_nse.geojson"),driver = driver)
    
    nodes_gdf_results = pd.merge(nodes_gdf, nse_nodes.reset_index(), left_on = 'compart_id', right_on = 'node').rename(columns={0 : 'nse'})
    nodes_gdf_results = pd.merge(nodes_gdf_results, infon_df.groupby(cluster).floodvol.sum().reset_index(), left_on = 'compart_id', right_on = cluster)
    nodes_gdf_results.to_file(os.path.join(results_root, "nodes_nse.geojson"),driver = driver)
