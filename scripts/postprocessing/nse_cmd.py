#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 17:15:14 2021

@author: bdobson
"""
import os
import pandas as pd
import geopandas as gpd
from glob import glob
import sys

if __name__ == "__main__":
    
    if len(sys.argv) > 1:
        runid = str(sys.argv[1])
        jobid = int(sys.argv[2])
        njobs = int(sys.argv[3])
    else:
        runid = '2021-03-28'
        jobid = 0
        njobs = 1

    root = os.path.join("/rds", "general", "user", "bdobson", "home", "cwsd_sewer","data")
    catchment = "cranbrook"
    sim_root = os.path.join(root,catchment,"results",runid)
    jobs = sorted(glob(os.path.join(sim_root,"*","*","*","*")))
    jobs = [x for x in jobs if os.path.isdir(x)]
    results_arcs = []
    results_nodes = []
    for (i, job) in enumerate(jobs):
        if (i % njobs) == jobid:
            path = os.path.normpath(job)
            job_split = path.split(os.sep)
            repetition = [x for x in job_split if 'repetition' in x][0]
            cluster_name = [x for x in job_split if 'cluster' in x][0]
            dt = [int(x.split('_')[2]) for x in job_split if 'dt' in x][0]
            rain = job_split[-1]
            N = int(cluster_name.split('_')[-1])
            
            
            compare_nodes = pd.read_csv(os.path.join(sim_root,repetition,"comparison_nodes.csv")).reset_index()[['node_id','quantile']]
            compare_arcs = pd.read_csv(os.path.join(sim_root,repetition,"comparison_arcs.csv")).reset_index()[['info_id','quantile']]
            
            node_gdf = gpd.read_file(os.path.join(sim_root,repetition,"clusters.geojson"))
            
            isjob = True
            flow_fid = os.path.join(job, "flows.gzip")
            if not os.path.isfile(flow_fid):
                isjob = False
            flood_fid = os.path.join(job, "flood_vol.gzip")
            if not os.path.isfile(flood_fid):
                isjob = False
            node_fid = os.path.join(job, "storages.gzip")
            if not os.path.isfile(node_fid):
                isjob = False
            
            if isjob:
                info_fid = os.path.join(job, "highfid_flows.gzip")
                
                info_df = pd.read_parquet(info_fid).set_index('time')
                info_df.index = pd.to_datetime(info_df.index)
                
                flow_df = pd.read_parquet(flow_fid)
                flow_df.index = pd.to_datetime(flow_df.index)
                
                info_fid = os.path.join(job, "highfid_nodes.gzip")
                
                infon_df = pd.read_parquet(info_fid).set_index('time')
                infon_df.index = pd.to_datetime(infon_df.index)
                
                node_df = pd.read_parquet(node_fid)
                node_df.index = pd.to_datetime(node_df.index)
                node_df.node = node_df.node.str.replace('-sewer','').astype(int)
                
                flood_df = pd.read_parquet(flood_fid)
                flood_df.index = pd.to_datetime(flood_df.index)
                flood_df.node = flood_df.node.str.replace('-land','').astype(int)
                
                flow_df['flow'] = flow_df['flow_out']
                
                mdf = pd.merge(node_df.reset_index(), flood_df.reset_index(),on=['date','node'],how='left').set_index('date')
                mdf = mdf.fillna(0)
                mdf['val'] = mdf.val_x + mdf.val_y
                
                flows = pd.merge(info_df[['arc','val']].reset_index(), flow_df[['arc','flow']].reset_index(), left_on = ['time','arc'], right_on = ['date','arc']).dropna()
                flows = flows.dropna(axis=0, how="any")
                # flows.loc[flows.val == 0, 'val'] = 0.001
                
                nodes = pd.merge(mdf[['node','val']].reset_index(), infon_df[[cluster_name,'volume']].reset_index(), left_on = ['date','node'], right_on = ['time', cluster_name]).dropna()
                # nodes.loc[nodes.val == 0, 'val'] = 0.001
                
                def f(x):
                    if x.val.var() > 0:
                        obj = 1 - sum((x.flow - x.val)**2)/sum((x.val - x.val.mean())**2)
                    else:
                        obj = None
                    return obj
                
                nse_flows = flows.loc[flows.arc.isin(compare_arcs.info_id)].groupby('arc').apply(f)
                
                def f(x):
                    if x.val.var() > 0:                    
                        obj = (x.flow - x.val).sum()/x.val.sum()
                    else:
                        obj = None
                    return obj
                
                bias_flows = flows.loc[flows.arc.isin(compare_arcs.info_id)].groupby('arc').apply(f)
                
                for arc in nse_flows.index:
                    results_arcs.append({'nse' : nse_flows.loc[arc],
                                         'bias' : bias_flows.loc[arc],
                                        'arc' : arc,
                                        'rain' : rain,
                                        'dt' : dt,
                                        'N' : N,
                                        'method' : cluster_name,
                                        'quantile' : compare_arcs.loc[compare_arcs.info_id == arc,'quantile'].values[0],
                                        'rep' : repetition
                                        })
                
                def f(x):
                    if x.val.var() > 0:
                        obj = 1 - sum((x.volume - x.val)**2)/sum((x.volume - x.volume.mean())**2)
                    else:
                        obj = None
                    return obj
                
                
                
                node_lookup = node_gdf[[cluster_name,'node_id']].drop_duplicates()
                node_lookup = node_lookup.loc[node_lookup.node_id.isin(compare_nodes.node_id)]
                nse_nodes = nodes.loc[nodes.node.isin(node_lookup[cluster_name])].groupby('node').apply(f)
                
                def f(x):
                    if x.val.var() > 0:                    
                        obj = (x.val - x.volume).sum()/x.volume.sum()
                    else:
                        obj = None
                    return obj
                bias_nodes = nodes.loc[nodes.node.isin(node_lookup[cluster_name])].dropna(how='any',axis=0).groupby('node').apply(f)
                # worst_floods = infon_df.dropna().groupby(cluster_name).max().volume
                
                for point in nse_nodes.index:
                    node_id = node_lookup.loc[node_lookup[cluster_name] == point,'node_id'].values[0]
                    results_nodes.append({'nse' : nse_nodes.loc[point],
                                          'bias' : bias_nodes.loc[point],
                                            'point' : point,
                                            'node_id' : node_id ,
                                            'rain' : rain,
                                            'dt' : dt,
                                            'N' : N,
                                            # 'maxstor' : worst_floods.loc[point],
                                            'method' : cluster_name,
                                            'quantile' : compare_nodes.loc[compare_nodes.node_id == node_id, 'quantile'].values[0],
                                            'rep' : repetition})
            
            
    results_nodes = pd.DataFrame(results_nodes)
    results_arcs = pd.DataFrame(results_arcs)
    
    fn = os.path.join(sim_root, "comparison_arcs_nse.csv")
    
    while True:	
        try:	
            os.rename(fn,fn)	
            results_arcs.to_csv(fn, mode = 'a',header=not os.path.exists(fn),index = False)	
            break	
        except:	
            pass	
    	
    fn = os.path.join(sim_root, "comparison_nodes_nse.csv")    	
    while True:	
        try:	
            os.rename(fn,fn)	
            results_nodes.to_csv(fn, mode = 'a',header=not os.path.exists(fn),index = False)	
            break	
        except:	
            pass
