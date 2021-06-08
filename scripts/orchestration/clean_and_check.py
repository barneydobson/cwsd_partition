#!/usr/bin/env python3

"""
Created on Tue Jun 16 15:19:24 2020

@author: Barney
"""

import os
import geopandas as gpd
import pandas as pd
from numpy import linspace, zeros, diff, where, NaN 
from shapely.geometry import Point, LineString
from shutil import copyfile, rmtree
import sys
from glob import glob
#Use to check all simulations have been completed
if __name__ == "__main__":
    
    if len(sys.argv) > 1:
        runid = str(sys.argv[1])
        jobid = int(sys.argv[2])
        njobs = int(sys.argv[3])
    else:
        runid = '2021-05-21'
        jobid = 0
        njobs = 1
    
    root = os.path.join("/rds","general","user","bdobson", "home", "cwsd_sewer","data")
    catchment = "cranbrook"
    sim_root = os.path.join(root,catchment,"results",runid)
    dts = [0.5,1,2,5,10,20,50,100,300,1440]
    # dts = [1000,1440]
    clusters = sorted(glob(os.path.join(sim_root,"*","*.geojson")))
 
    files = []
    roots = []
    print('generating job scheduler')
    for cluster in clusters:
          cols = gpd.read_file(cluster).columns
          
          cols = cols[cols.str.contains('cluster')]
          for col in cols:
              if "L_L_sc" in col:
                  dts_ = dts
              else:
                  dts_ = [1]
              for dt in dts_:
                  fid = os.path.join(os.path.split(cluster)[0],
                                         col,
                                         '_'.join(['sim_dt',str(int(dt*60)),'s']))
                  roots.append(fid)
                  for rain in ['august','january']:
                    for fn in ['flows.gzip','depths.gzip','flood_vol.gzip','storages.gzip']:
                      files.append(os.path.join(fid,rain,fn))
    print(roots[0])
    print(roots[-1])
    print('reading run jobs')                    
    run_jobs = glob(os.path.join(sim_root,"*","*","*"))
    run_files = []
    for run_job in run_jobs:
        if run_job not in roots:
          #rmtree(run_job)
          print(run_job)
        else:
          rains = [ f.name for f in os.scandir(run_job) if f.is_dir() ]
          for rain in rains:
            for fn in os.scandir(os.path.join(run_job, rain)):
              run_files.append(os.path.join(run_job,rain,fn.name))
    print('diff')
    missing = []
    for item in files:
      if item not in run_files:
        missing.append(os.path.split(os.path.split(item)[0])[0]) 
    pd.Series(missing).drop_duplicates().to_csv('missing_sims_new.csv')

