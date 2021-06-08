# -*- coding: utf-8 -*-
"""
Created on Thu Jan 14 13:41:46 2021

@author: bdobson
"""

import os
import pandas as pd
import geopandas as gpd

PI = 3.14159265359
driver = 'GeoJSON'
extension = '.geojson'
SMALL_NUMBER = 1e-12 # Small but detectable number to prevent 0 divide by 0
LARGE_NUMBER = 1e12 # Large number for inf check
#Addresses
data_root = os.path.join("C:\\", "Users", "bdobson", "Documents", "GitHub", "cwsd_sewer", "data", "cranbrook")

processed_root = os.path.join(data_root, "processed")
raw_root = os.path.join(data_root, "raw")

info_dir = "info_sims"
back_calc_storm = "january"

#Load data
edges_gdf = gpd.read_file(os.path.join(processed_root, "edges_clean.geojson"))
nodes_gdf = gpd.read_file(os.path.join(processed_root, "nodes_clean.geojson"))

flow_df = pd.read_csv(os.path.join(raw_root,info_dir,back_calc_storm, "ds_flow.csv")).set_index(['Time','Seconds'])
flow_df = flow_df.reindex(columns=edges_gdf.info_id)

vel_df = pd.read_csv(os.path.join(raw_root,info_dir,back_calc_storm, "ds_vel.csv")).set_index(['Time','Seconds'])
vel_df = vel_df.reindex(columns=edges_gdf.info_id)

us_vel_df = pd.read_csv(os.path.join(raw_root,info_dir,back_calc_storm, "us_vel.csv")).set_index(['Time','Seconds'])
us_vel_df = us_vel_df.reindex(columns=edges_gdf.info_id)

dep_df = pd.read_csv(os.path.join(raw_root,info_dir,back_calc_storm, "depnod.csv")).set_index(['Time','Seconds'])
stor_df = pd.read_csv(os.path.join(raw_root,info_dir,back_calc_storm, "volume.csv")).set_index(['Time','Seconds'])
flood_df = pd.read_csv(os.path.join(raw_root,info_dir,back_calc_storm, "floodvolume.csv")).set_index(['Time','Seconds'])

ind = (flow_df > 0) & (vel_df <= 0)
vel_df[ind] = us_vel_df[ind]

vel_df['node_1656.2'] = 0.75 # m/s (infoworks manual - can be between 0.75-1.8 for pumps)

xc_df = flow_df.div(vel_df).max()
xc_df[(xc_df < 0) | xc_df.isna() | (xc_df > LARGE_NUMBER) ] = 0
xc_df = xc_df.reset_index().rename(columns={'index' : 'info_id', 0 : 'cross_sec'})
flow_df = flow_df.max().reset_index().rename(columns={'index' : 'info_id', 0 : 'capacity'})

#Back calculate missing XC and capacity from simulation data
edges_gdf.loc[edges_gdf.cross_sec.isna(),'cross_sec'] = xc_df.loc[edges_gdf.cross_sec.isna()]

edges_gdf.loc[edges_gdf.capacity.isna(),'capacity'] = flow_df.loc[edges_gdf.capacity.isna()]
edges_gdf.loc[edges_gdf.capacity <= 0,'capacity'] = flow_df.loc[edges_gdf.capacity <= 0]

#Back calculate missing length from geometry data
edges_gdf.loc[edges_gdf['length'].isna(), 'length'] = edges_gdf.geometry.length.loc[edges_gdf['length'].isna()]

#Mannings
edges_gdf['mannings'] = 70 # (the rivers are also 70)

#Back calculate missing slope assuming that everything is a circular pipe in full bore
r = (edges_gdf.cross_sec / PI) ** 0.5
R = edges_gdf.cross_sec / (2 * PI * r)

slope = (edges_gdf.capacity / (edges_gdf.cross_sec * (R ** (2/3)) * edges_gdf.mannings)) ** 2
slope[(slope < 0) | slope.isna() | (slope > LARGE_NUMBER)] = 0

edges_gdf.loc[edges_gdf.gradient <= 0,'gradient'] = slope.loc[edges_gdf.gradient <= 0]
edges_gdf.loc[edges_gdf.gradient.isna(),'gradient'] = slope.loc[edges_gdf.gradient.isna()]

#Back calculate node storages for non-sensible storages
nodes_gdf = pd.merge(nodes_gdf,(stor_df - flood_df).max().rename('storage_sim').reset_index(), left_on = 'node_id', right_on = 'index',how='left')

nodes_gdf.loc[nodes_gdf.storage <= 0, 'storage'] = nodes_gdf.loc[nodes_gdf.storage <= 0, 'storage_sim']
nodes_gdf = nodes_gdf.drop(['index','storage_sim'],axis=1)

nodes_gdf.loc[nodes_gdf.chamber_ar == 0,'chamber_ar'] = 0.5

#Really odd weirs
node_ind = 'node_2477'
nodes_gdf.loc[nodes_gdf.node_id == node_ind,'storage'] = (dep_df[node_ind].max() - nodes_gdf.loc[nodes_gdf.node_id == node_ind,'chamber_fl']) * nodes_gdf.loc[nodes_gdf.node_id == node_ind,'chamber_ar']

node_ind = 'lake_3'
nodes_gdf.loc[nodes_gdf.node_id == node_ind,'chamber_fl'] = dep_df[node_ind].min()
nodes_gdf.loc[nodes_gdf.node_id == node_ind,'chamber_ar'] = (stor_df[node_ind]/(dep_df[node_ind] - (dep_df[node_ind].min() - 0.01))).median()

node_ind = 'lake_2'
nodes_gdf.loc[nodes_gdf.node_id == node_ind,'chamber_fl'] = dep_df[node_ind].min()

edges_gdf.to_file(driver = driver, filename= os.path.join(processed_root,"edges_clean" + extension))
nodes_gdf.to_file(driver = driver, filename= os.path.join(processed_root,"nodes_clean" + extension))
