# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 14:37:08 2020

@author: bdobson
"""
#warning... this script is a horrible mess because no infoworks model is perfect

import os
import geopandas as gpd
import pandas as pd
import numpy as np
from shapely.geometry import LineString, Point

def clean(data_root:str,
          catchment='cranbrook' ):
    CRS = "EPSG:27700"
    HA_TO_M2 = 10000
    PCT_TO_PROP = 0.01
    MM_TO_M = 0.001
    PI = 3.14159265359
    G = 9.81
    # catchment = "cranbrook"
    driver = 'GeoJSON'
    extension = '.geojson'
    
    raw_root = os.path.join(data_root,"raw")
    output_root = os.path.join(data_root,"processed")
    xl_data = os.path.join(raw_root, "Network data.xlsx")
    
    #Load data
    def load_arc(title):
        arcs_fid = os.path.join(raw_root,title)
        arcs_gdf = gpd.read_file(arcs_fid)
        arcs_gdf[['us_node_id','ds_node_id']] = arcs_gdf[['us_node_id','ds_node_id']].astype(str)
        arcs_gdf['info_id'] = arcs_gdf.us_node_id + '.' + arcs_gdf.link_suffi
        arcs_gdf['edge_type'] = title.replace('.shp','').lower()[0:-1]
        if title == 'Weirs.shp':

                arcs_gdf['weir_height'] = arcs_gdf['crest']

                #Assumes full upstream and zero downstream level
                arcs_gdf['weir_param_mul'] = (G ** 0.5) * \
                                                 arcs_gdf['discharge_'] * \
                                                 arcs_gdf['width']
                arcs_gdf['weir_param_pow'] = 1.5                                 
                
                # theta = arcs_gdf['notch_angl'] * PI / 180
                # arcs_gdf.loc[~ind,'capacity'] = (8/15) *\
                #                                   (h.loc[~ind] ** (5/2)) * \
                #                                   np.tan(theta/2) * \
                #                                   ((2 * G) ** 0.5) * \
                #                                   arcs_gdf.loc[~ind,'discharge_']

        elif title == 'Orifices.shp':
            a = arcs_gdf['diameter']**2 * PI /4
            arcs_gdf['capacity'] = arcs_gdf['discharge_'] * \
                                     a * \
                                     (G ** 0.5) * \
                                         (arcs_gdf['invert'] ** 0.5)
        
        return arcs_gdf
    links_gdf = load_arc('Conduits.shp')
    weirs_gdf = load_arc('Weirs.shp')
    orifices_gdf = load_arc('Orifices.shp')
    pumps_gdf = load_arc('Pump.shp')
    rivers_gdf = load_arc('River_reaches.shp')
    links_gdf = pd.concat([links_gdf, weirs_gdf, orifices_gdf, pumps_gdf, rivers_gdf])
    
    
    def load_node(title):
        nodes_fid = os.path.join(raw_root,title)
        nodes_gdf = gpd.read_file(nodes_fid)
        nodes_gdf['node_id'] = nodes_gdf['node_id'].astype(str)
        if title == 'Sub-catchments_polygons.shp':
            nodes_gdf = nodes_gdf.rename(columns = {'catchment_' : 'slope'})
        return nodes_gdf
    
    nodes_gdf = load_node("Nodes.shp")
    catch_gdf = load_node("Sub-catchments_polygons.shp")
    nodes_gdf = pd.merge(nodes_gdf, catch_gdf, on = 'node_id', how = 'outer')
    nodes_gdf = nodes_gdf.rename(columns={'geometry_x' : 'geometry', 
                                          'geometry_y' : 'subcatch_geom',
                                          'area_perce' : 'area_perc0',
                                          'area_per01' : 'area_perc10',
                                          'area_per02' : 'area_perc11'})
    del catch_gdf, weirs_gdf, pumps_gdf, orifices_gdf, rivers_gdf
    
    
    if catchment == 'cranbrook':
        #Suspected errors
        nodes_gdf.loc[nodes_gdf.node_id == "node_2098", 'node_type'] = "Manhole"
        nodes_gdf = nodes_gdf.loc[~nodes_gdf.node_id.isin(["node_412","node_418"])]
 
    #Here we assume that if there are 2 sub-catchments for the same node_id they can be aggregated together
    def f(x):
        d = {}
        d['total_area'] = x['total_area'].sum()
        d['slope'] = (x['total_area'] * x['slope'] / d['total_area']).sum()
        for y in range(12):
            ind = 'area_perc' + str(y)
            d[ind] = (x['total_area'] * x[ind] / d['total_area']).sum()
        d['chamber_fl'] = x['chamber_fl'].mean()
        d['chamber_ar'] = x[['chamber_ar','base_area']].max(axis=1).sum()
        d['chamber_ro'] = x['chamber_ro'].mean()
        d['flood_leve'] = x['flood_leve'].mean()
        d['node_type'] = x['node_type'].iloc[0]
        d['geometry'] = Point((x.geometry.x).mean(),(x.geometry.y).mean())
        
        return pd.Series(d)
    
    nodes_gdf = gpd.GeoDataFrame(nodes_gdf.groupby('node_id').apply(f).reset_index())

        
    #Calculate runoff info for nodes
    nodes_clean = nodes_gdf[['node_id','node_type','geometry','chamber_fl','chamber_ar']].copy()

    if catchment == 'cranbrook':
        run_sheetname = 'Runoff surface'
        runoff_surfs = pd.read_excel(xl_data, sheet_name = run_sheetname).set_index('Runoff surface ID')
    
    nodes_clean['storage'] = (nodes_gdf['flood_leve'].fillna(nodes_gdf['chamber_ro']) - nodes_gdf['chamber_fl'])*nodes_gdf['chamber_ar']
    # nodes_clean['floor_ad'] = nodes_gdf['chamber_fl']
    nodes_clean['area'] = nodes_gdf['total_area'] * HA_TO_M2
    nodes_clean['run_coef'] = 0
    nodes_clean['init_loss'] = 0
    
    if catchment == 'cranbrook':
        lakes = {'lake_0' : 967.5,
                    'lake_1' : 4407.63,
                    'lake_2' : 23482,
                    'lake_3' : 48676,
                    'lake_4' : 140.82}
        nodes_clean = nodes_clean.set_index('node_id')
        nodes_clean.loc[lakes.keys(),'storage'] = pd.Series(lakes)
        nodes_clean = nodes_clean.reset_index()
    
        for surface in nodes_gdf.loc[:,nodes_gdf.columns.str.contains('area_perc')].columns:
            surface_id = int(surface.replace('area_perc','')) + 1
            coef = runoff_surfs.loc[surface_id, 'Fixed runoff coefficient']
            loss = runoff_surfs.loc[surface_id, 'Initial loss value (m)']
            area = nodes_gdf[surface]
            nodes_clean['run_coef'] += coef * area * PCT_TO_PROP
            nodes_clean['init_loss'] += loss * area * PCT_TO_PROP
        
    #Remove any null geometries
    nodes_clean = nodes_clean.loc[nodes_clean.geometry.x.notna()]
    
    #Calculate volume/flow info for links
    if catchment == 'cranbrook':
        edges_clean = links_gdf[['us_node_id','ds_node_id','conduit_le','capacity','gradient','geometry','info_id','edge_type','weir_param_pow','weir_param_mul','weir_height']].copy().rename(columns={'conduit_le' : 'length'})
    
    edges_clean['mannings'] = links_gdf[["bottom_ro1","top_rough1"]].mean(axis=1)
    ind_circ = links_gdf.shape_1 == 'CIRC'
    edges_clean.loc[ind_circ,'cross_sec'] = links_gdf.loc[ind_circ, 'conduit_wi'].mul(MM_TO_M).pow(2)*PI/4
        
    edges_clean.loc[~ind_circ,'cross_sec'] = links_gdf.loc[~ind_circ, 'conduit_wi'].mul(MM_TO_M) * links_gdf.loc[~ind_circ, 'conduit_he'].mul(MM_TO_M)
    edges_clean = edges_clean.loc[edges_clean.geometry.notna()]
    
    #Double check why I need to do this
    edges_clean.capacity = edges_clean.capacity.astype(float)
    
    #Write cleaned dataframes
    edges_clean.crs = CRS
    nodes_clean.crs = CRS
    nodes_clean.to_file(driver = driver, filename= os.path.join(output_root,"nodes_clean" + extension))
    edges_clean.to_file(driver = driver, filename= os.path.join(output_root,"edges_clean" + extension))

catchment = "cranbrook"
data_root = os.path.join("C:\\", "Users", "bdobson", "Documents", "GitHub", "cwsd_sewer", "data", "cranbrook")

clean(data_root)