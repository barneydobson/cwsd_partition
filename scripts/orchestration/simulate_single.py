#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 22 10:58:27 2021

@author: bdobson
"""

import os
import geopandas as gpd
import pandas as pd
import numpy as np
import models
import constants
from glob import glob
from shutil import copyfile
from datetime import datetime
import sys
def sim(cluster_root:str,
        rain_fid:str,
        surrogate_type:str='physical',
        LAG:int=200,
        dt_sim = 5):
    """ Runs the model.
        Parameters: surrogate_type can be 'physical', 'hybrid', 'cwsd'
                    LAG should be 102 for 30 mins rain or 0 for 960 mins
    """
    """Additional constants
    """
    DT_INTENSITY = pd.Timedelta('1h') #i.e. assumes rainfall is in per hour intensity
    if '0min' in rain_fid:
        DT_INPUT = pd.Timedelta('1m') #Needs to match the same DT in inputs
    else:
        DT_INPUT = pd.Timedelta('5m') #Needs to match the same DT in inputs
        
    DT_SIM = pd.Timedelta(value = dt_sim, unit = 'm') #Needs to match the same DT_SIM in orchestration

    LAG = 200
    constants.MM_M2_TO_SIM_VOLUME = 1e-3 # m3
    constants.M3_S_TO_M3_DT = DT_SIM / pd.Timedelta('1s')
    
    """Options
    """
    # catchment = "cranbrook" # can be 'cranbrook'
    # rain_name = "rain_960min"
    extension = '.geojson'
    sim_start = '2000-01-01' #Misc start date that detailed simulation results are assumed to start at if actual date isn't given
    
    """Load and convert data
    """
    nodes_fid = os.path.join(cluster_root,"compartment_nodes" + extension)
    nodes_gdf = gpd.read_file(nodes_fid).apply(pd.to_numeric,errors='ignore').set_index('compart_id')
    nodes_gdf.index = nodes_gdf.index.astype(int)
    nodes_gdf = nodes_gdf.sort_values(by='upstr_ind', ascending = True).rename(columns={'storage' : 'capacity'})
    
    nodes_gdf.loc[nodes_gdf.chamber_ar == 0, "chamber_ar"] = 0.5
    
    edges_fid = os.path.join(cluster_root,"compartment_edges" + extension)
    edges_gdf = gpd.read_file(edges_fid).apply(pd.to_numeric,errors='ignore')
    
    nodes_gdf['compartment_type'] = 'basic'
    nodes_gdf.loc[~nodes_gdf.index.isin(edges_gdf.start_comp),'compartment_type'] = 'outfall'
    
    
    storm_df = pd.read_csv(rain_fid, sep=',')
    padding = pd.DataFrame(index = range(len(storm_df),len(storm_df)+LAG), data=np.zeros(LAG), columns = ['intensity_mm_hr'])
    storm_df = storm_df.append(padding)
    
    if 'date' in storm_df.columns:
        sim_start = pd.to_datetime(storm_df.iloc[0].date)
    
    storm_df['rain'] = storm_df.intensity_mm_hr * DT_INPUT / DT_INTENSITY
    
    storm_df.index = pd.date_range(start = pd.to_datetime(sim_start), periods = len(storm_df), freq = DT_INPUT)
    if DT_SIM >= DT_INPUT:
        storm_df = storm_df.rain.resample(DT_SIM).sum()
    else:
        storm_df = storm_df.rain.resample(DT_SIM).interpolate()*(DT_SIM/DT_INPUT)
        
    storm_df.index = storm_df.index.astype(str)
    storm_df = storm_df.to_dict()
        
    my_model = models.model(physical = True)
    my_model.add_nodes(nodes_gdf.drop('geometry',axis=1).T.to_dict())
    my_model.add_arcs(edges_gdf.drop('geometry',axis=1).set_index('info_id').T.to_dict())
    my_model.add_inputs({'rain_holder' : storm_df})  
    
    

    for node in my_model.model_nodes.values():
        node.model = my_model
    for arc in my_model.model_arcs.values():
        arc.model = my_model
    
    results = my_model.run()
    return results
        
if __name__ == "__main__":
    
   
    runid = '2021-06-07'
    jobid = 0
    njobs = 1

    root = os.path.join("C:\\", "Users", "bdobson", "Documents", "GitHub", "cwsd_sewer","data")
    catchment = "cranbrook"
    clusters_root = os.path.join(root,catchment,"results",runid)
    
    times = []
    jobs = glob(os.path.join(clusters_root,"*","*"))
    
    for (i, cluster) in enumerate(jobs):
        dt = int(os.path.split(cluster)[-1].split('_')[-2]) / 60
        rains = [ f.name for f in os.scandir(cluster) if f.is_dir() ]
        
        cluster_type = os.path.split(os.path.split(cluster)[0])[-1]
        rep = os.path.split(os.path.split(cluster)[0])[0]
        
        for rain in rains:
            print(cluster,rain)
            fid = os.path.join(root, catchment, "raw", "storms", rain + '.csv')
            start = datetime.now()
            results = sim(cluster, fid, dt_sim = dt)
            end = datetime.now()
            times.append({'fid' : cluster,
                          'N' : int(os.path.split(os.path.split(cluster)[0])[-1].split('_')[-1]),
                          'rain' : rain,
                          'dt' : dt,
                          'time_s' : (end - start).total_seconds()})
            
        
            for key, value in results.items():
                value = pd.DataFrame(value)
                if key == 'flows':
                    value.set_index('date').to_parquet(os.path.join(cluster,rain, key + '.gzip'), compression = 'GZIP')
                else:
                    nodes = value.node.str.replace('-land','').str.replace('-sewer','').astype(int)
                    value.set_index('date').to_parquet(os.path.join(cluster,rain, key + '.gzip'), compression = 'GZIP')
            
            del results
    
            for fid in [__file__, "models.py"]:
                fn = os.path.basename(fid)
            copyfile(fid, os.path.join(cluster, rain, fn))
# times = pd.DataFrame(times)
# times.to_csv(os.path.join(clusters_root, "runtimes.csv"), index = False)

