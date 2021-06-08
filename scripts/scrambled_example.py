# -*- coding: utf-8 -*-
"""
Created on Tue Jun  8 17:01:52 2021

@author: bdobson
"""

import os
import sys
import geopandas as gpd
import pandas as pd
import networkx as nx
import numpy as np

repo_root = os.path.join("C:\\", "Users", "bdobson", "Documents", "GitHub", "cwsd_partition")
sys.path.append(os.path.join(repo_root, "scripts","preprocessing"))

import partitioning
import aggregate_cmd

#Load data
edges_gdf = gpd.read_file(os.path.join(repo_root, "data", "cranbrook", "processed", "edges_gdf_scrambled.geojson"))
nodes_gdf = gpd.read_file(os.path.join(repo_root, "data", "cranbrook", "processed", "nodes_gdf_scrambled.geojson"))
coords = pd.DataFrame(index = nodes_gdf.node_id, data = np.array( [[i.x, i.y] for i in nodes_gdf.geometry] ))


#Create networkx graph
edges = [(x.us_node_id, x.ds_node_id) for x in edges_gdf.itertuples()]

G = nx.DiGraph()
G.add_nodes_from(nodes_gdf.node_id)
G.add_edges_from(edges)
del edges


#Perform clustering
N_clusters = 100 # Will be ignored if using 'Louv'
method = 'louv'
x = partitioning.wrapper(nx.to_scipy_sparse_matrix(G), 0, coords, method)
nodes_gdf[method] = x


