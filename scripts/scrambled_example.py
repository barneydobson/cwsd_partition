# -*- coding: utf-8 -*-
"""
Created on Tue Jun  8 17:01:52 2021

@author: bdobson
"""

import os
import sys
import geopandas as gpd
import pandas as pd
import datetime
import networkx as nx
import numpy as np

#Load CWSD modules
repo_root = os.path.join("C:\\", "Users", "Barney", "Documents", "GitHub", "cwsd_partition")
sys.path.append(os.path.join(repo_root, "scripts","preprocessing"))

import partitioning
import aggregate_cmd

sys.path.append(os.path.join(repo_root, "scripts","orchestration"))
import simulate_cmd

#Specify rainfall timeseries (see data/cranbrook/raw)
rain_fid = os.path.join(repo_root, "data", "cranbrook", "raw", "storms", "august.csv")

#Load data
edges_gdf = gpd.read_file(os.path.join(repo_root, "data", "cranbrook", "processed", "edges_gdf_scrambled.geojson"))
nodes_gdf = gpd.read_file(os.path.join(repo_root, "data", "cranbrook", "processed", "nodes_gdf_scrambled.geojson"))
coords = pd.DataFrame(index = nodes_gdf.node_id, data = np.array( [[i.x, i.y] for i in nodes_gdf.geometry] ))

#Make results directory structure
results_name = datetime.datetime.today().date().strftime(format='%Y-%m-%d')
results_folder = os.path.join(repo_root, "data", "cranbrook", "results", results_name)
if not os.path.exists(results_folder):
    os.mkdir(results_folder)

#Create networkx graph
edges = [(x.us_node_id, x.ds_node_id) for x in edges_gdf.itertuples()]

G = nx.DiGraph()
G.add_nodes_from(nodes_gdf.node_id)
G.add_edges_from(edges)
del edges

#Perform partitioning
N_partitions = 100 # Will be ignored if using 'louv'
method = 'louv'
x = partitioning.wrapper(nx.to_scipy_sparse_matrix(G), N_partitions, coords, method)
nodes_gdf[method] = x

#Preserve outfalls by putting them in their own partition (aggregation assumes this)
outfalls = list(set(nodes_gdf.node_id).difference(edges_gdf.us_node_id))
nodes_gdf.loc[nodes_gdf.node_id.isin(outfalls), method] = np.arange(nodes_gdf[method].max() + 1, nodes_gdf[method].max() + len(outfalls) + 1)

#Save partitions
nodes_gdf.to_file(driver = 'GeoJSON', filename = os.path.join(results_folder, 'partitions.geojson'))

#Specify timestep
dt = 1 # (minutes)

#Aggregate to physical model (you can view the reduced model in results_folder)
aggregate_cmd.aggregate(os.path.join(repo_root, "data"),
                        results_folder,
                        method,
                        dt_sim = dt, 
                        demo = True)

#Simulate
results = simulate_cmd.sim(os.path.join(results_folder,
                                        method,
                                        "_".join(["sim","dt",str(int(dt*60)),"s"])),
                           rain_fid,
                           dt_sim = dt)