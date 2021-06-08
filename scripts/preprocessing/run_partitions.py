# -*- coding: utf-8 -*-
"""
Created on Fri Jul 24 17:56:13 2020

@author: hanna
"""
import os
import geopandas as gpd
import pandas as pd
import numpy as np
import networkx as nx
import datetime

from shutil import copyfile
import partitioning

#Calls different partitioning algorithms

start_all = datetime.datetime.now()
def create_clusters(desired_n: list,
                    methods: list,
                    root: str,
                    root_results: str,
                    location: str = 'cranbrook',
                    weighting: str = 'capacity',
                    repetition: int = 0
                    ):
    """ Performs 7 different clustering algorithms on the data set for each n.
        methods can contain: 'sc','ac','km','ac_sc','km_sc','L_scn','L_L_sc','louv'
        Returns geojsons of nodes and clusters for each n and each algorithm,
        as well as a cvs file of metrics for each algorithm.
    """
    
    
    repetition_results = os.path.join(root_results, "repetition_" + str(repetition))
    if not os.path.exists(repetition_results):
        os.mkdir(repetition_results)
    """Load
    """
    edges_fid = os.path.join(root, "processed", "edges_clean.geojson")
    edges_gdf = gpd.read_file(edges_fid)

    nodes_fid = os.path.join(root, "processed", "nodes_clean.geojson")
    nodes_gdf = gpd.read_file(nodes_fid)
    results_gdf = nodes_gdf.set_index('node_id').copy()
    
    coords = pd.DataFrame(index = nodes_gdf.node_id, data = np.array( [[i.x, i.y] for i in nodes_gdf.geometry] ))
    
    info_flow_fid = os.path.join(root, "raw", "info_sims", "january", "ds_flow.csv")
    info_node_fid = os.path.join(root, "raw", "info_sims", "january", "floodvolume.csv")
    
    
    def sample_info(fid):
        df = pd.read_csv(fid)
        df = df.drop(['Time','Seconds'],axis=1).sum()
        df = df.loc[df > 0]
        preserve_quantiles = [1]
        samples_per_quantile = 1
        
        preserved = [df.loc[(df < x) & (df > y)].sample(samples_per_quantile) for x, y in zip(df.quantile(preserve_quantiles), df.quantile([0.8] + preserve_quantiles[:-1]))]
        preserved = pd.concat(preserved)
        preserved.loc[:] = np.repeat(preserve_quantiles, samples_per_quantile)
        # preserved = {x : y.index.tolist() for x, y in zip(preserve_quantiles, preserved)}
        return preserved
    
    comparison_nodes = sample_info(info_node_fid)
    comparison_arcs = sample_info(info_flow_fid)
    
    #Save comparison points
    comparison_nodes.rename('quantile').reset_index().rename(columns={'index':'node_id'}).to_csv(os.path.join(repetition_results, "comparison_nodes.csv"))
    comparison_arcs.rename('quantile').reset_index().rename(columns={'index':'info_id'}).to_csv(os.path.join(repetition_results, "comparison_arcs.csv"))
    del nodes_fid, edges_fid
    
    
    """Force clusters
    """
    #Keep lakes, outfalls and starting nodes for orifices and weirs out of clusters

    # ind = nodes_gdf.node_type.isin(['Outfall']) 
    ind = nodes_gdf.node_type.isin(['Outfall','Storage']) 
    ind = ind | (nodes_gdf.node_id.isin(edges_gdf.loc[edges_gdf.edge_type.isin(['weir']),'us_node_id']))
    ind = ind | (nodes_gdf.node_id.isin(edges_gdf.loc[edges_gdf.edge_type.isin(['weir']),'ds_node_id']))
    
    #Preserve comparison points
    ind = ind | nodes_gdf.node_id.isin(comparison_nodes.index)
    ind = ind | (nodes_gdf.node_id.isin(edges_gdf.loc[edges_gdf.info_id.isin(comparison_arcs.index),'us_node_id']))
    ind = ind | (nodes_gdf.node_id.isin(edges_gdf.loc[edges_gdf.info_id.isin(comparison_arcs.index),'ds_node_id']))
    
    preserved_nodes = nodes_gdf.loc[ind]
    nodes_gdf = nodes_gdf.loc[~ind]
    edges_gdf = edges_gdf.loc[edges_gdf.us_node_id.isin(nodes_gdf.node_id) & edges_gdf.ds_node_id.isin(nodes_gdf.node_id)]
    del ind
    
    
    """Create Graph and subgraphs
    """
    edges = [(x.us_node_id, x.ds_node_id) + ((getattr(x, weighting),) if weighting else tuple()) for x in edges_gdf.itertuples()]

    G = nx.DiGraph()
    G.add_nodes_from(nodes_gdf.node_id)
    G.add_weighted_edges_from(edges) if weighting else G.add_edges_from(edges)
    del edges
    G_ = []
    for (i,sg) in enumerate(nx.weakly_connected_components(G)):
        sg = G.subgraph(sg).copy()
        G_.append({'graph_ind' : i,
                   'nxgraph' : sg,
                   'n_nodes' : len(sg.nodes),
                   'adj_mat' : nx.to_scipy_sparse_matrix(sg),
                   'coords'  : coords.loc[list(sg.nodes)].values
                   })
        
        
    tot_nodes = len(G.nodes)
    
    """Run clustering
    """
    times = []
    for method in methods:
        for n in desired_n:
            print(method,n)
            start=datetime.datetime.now()
            sg_results = []
            for sg in G_:
                nn = int(n * sg['n_nodes'] / tot_nodes)
                if nn >= sg['n_nodes']:
                    x = np.arange(sg['n_nodes'])
                elif nn == 0:
                    x = [0] * sg['n_nodes']
                else:
                    x = partitioning.wrapper(sg['adj_mat'], nn, sg['coords'], method)
                sg_results.append({'nodes' : list(sg['nxgraph'].nodes), 'partitions' : x})
                
            end=datetime.datetime.now()
            times.append({'method' : method,
                          'n' : n,
                          'ns' : len(G_),
                          'time' : end - start})
            for nid in preserved_nodes.node_id:
                sg_results.append({'nodes' : [nid], 'partitions' : [0]})
            
            g_nodes = []
            g_partitions = []
            num_p = 0
            for sg in sg_results:
                g_nodes += sg['nodes']
                g_partitions += [x + num_p for x in sg['partitions']]
                num_p = max(g_partitions) + 1
            namestr = '_'.join(['cluster', method,'n',str(len(set(g_partitions)))])
            
            results_gdf[namestr] = pd.Series(index = g_nodes, data = g_partitions, name = namestr)
    
    
    times = pd.DataFrame(times)
    times.time = times.time.dt.total_seconds()
    results_gdf.to_file(driver = 'GeoJSON', filename = os.path.join(repetition_results, 'clusters.geojson'))
    times.to_csv(os.path.join(repetition_results,"cluster_timings.csv"),index=False)
    
    
    return results_gdf

CRS = "EPSG:27700"
location = 'cranbrook'

root = os.path.join("C:\\", "Users", "bdobson", "Documents", "GitHub", "cwsd_sewer", "data", location)

results_name = datetime.datetime.today().date().strftime(format='%Y-%m-%d')
root_results = os.path.join(root, "results", results_name)
if not os.path.exists(root_results):
    os.mkdir(root_results)
np.random.seed(10)

fid = __file__
fn = os.path.basename(fid)
copyfile(fid, os.path.join(root_results, fn))

methods = ['louv', "L_L_sc"]
# methods = ["L_L_sc","km_sc","L_scn","fluid","ac_sc","louv","sc","ac","km",]
desired_n = [500,1,5,10,20,50,100,200]
# desired_n = [200]
weighting = 'capacity'
repetitions = 20
for repetition in range(repetitions):
    
    print('rep = ' + str(repetition))
    create_clusters(desired_n, methods, root, root_results, location, weighting,repetition)