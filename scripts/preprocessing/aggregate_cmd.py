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
from shutil import copyfile
import sys
from glob import glob

def aggregate(root: str,
              cluster_root: str,
              cluster: str,
              catchment: str = "cranbrook",
              dt_sim = 5):
    """ Aggregates clusters into compartments.
        Paremeters: catchment can be "cranbrook" or "norwich"
        Returns a file 
    """
    E = 2.71828
    
    DT_ARC = pd.Timedelta('1s') #i.e. assumes arc flows are given in average m3/s over timestep
    DT_RAW = pd.Timedelta('1m') #Timestep of the infoworks files (both shp and results)
    DT_SIM = pd.Timedelta(value = dt_sim, unit = 'm') #Needs to match the same DT_SIM in orchestration
    M3_S_TO_M3_DT = DT_SIM / pd.Timedelta('1s') 
    S_TO_DT = 1 / DT_SIM.total_seconds() #Convert seconds to timestep
    SURFACE_TRAVEL_TIME = 5 * 60 # seconds
    TIMEAREA_RESOLUTION = 1000 #Resolution (s) at which timearea diagrams of individual nodes are superimposed at
    CRS = "EPSG:27700"
    # catchment = "cranbrook" # can be 'cranbrook' or 'norwich'
    driver = 'GeoJSON'
    extension = '.geojson'
    sim_start = '2000-01-01' #Misc start date that detailed simulation results are assumed to start at if actual date isn't given
    
    sim_rains = ["january","august"]
    
    """Load data
    """
    data_root = os.path.join(root, catchment,"processed")
    raw_root = os.path.join(root, catchment,"raw")
    
    edges_fid = os.path.join(data_root,"edges_clean.geojson")
    edges_gdf = gpd.read_file(edges_fid)
    edges_gdf.capacity = edges_gdf.capacity.astype(float)
    
    nodes_fid = os.path.join(cluster_root,"clusters.geojson")
    nodes_gdf = gpd.read_file(nodes_fid).rename(columns={'area' : 'surf_area'})

    info_dir = "info_sims"
    info_dir = os.path.join(raw_root,info_dir)    
    
    #InfoWorks node simulation for initial storage only    
    rain = "january"
    info_volume_fid = os.path.join(info_dir,rain, "volume.csv")
    info_volume_df = pd.read_csv(info_volume_fid)
    del nodes_fid, edges_fid, info_volume_fid
        
    edges_gdf['travel_time'] = 0
    ind = edges_gdf.capacity > 0
    edges_gdf.loc[ind,'travel_time'] = edges_gdf.loc[ind,'length'] * edges_gdf.loc[ind, 'cross_sec'] / edges_gdf.loc[ind,'capacity'] #s
    edges_gdf.capacity *= M3_S_TO_M3_DT
    
    """Manual corrections
    """   
    nodes_gdf = nodes_gdf.set_index('node_id').join(info_volume_df.iloc[0]).rename(columns={info_volume_df.iloc[0].name : 'initial_storage'}).reset_index()
    nodes_gdf.initial_storage = nodes_gdf.initial_storage.fillna(0)
    # nodes_gdf.storage = nodes_gdf.storage - nodes_gdf.initial_storage
    # nodes_gdf = nodes_gdf.drop('initial_storage', axis=1)
    
    """Classify edges within compartments
    """
    edges_gdf['storage'] = edges_gdf['length'] * edges_gdf['cross_sec']

    
    for cluster_type in [cluster]:
        
        cluster_folder = os.path.join(cluster_root, cluster_type)
        if not os.path.exists(cluster_folder):
            os.mkdir(cluster_folder)
        cluster_folder = os.path.join(cluster_folder, "_".join(["sim","dt",str(int(DT_SIM.total_seconds())),"s"]))
        if not os.path.exists(cluster_folder):
            os.mkdir(cluster_folder)
        
        
        edges_cluster = pd.merge(edges_gdf, nodes_gdf[['node_id',cluster_type]],left_on='us_node_id',right_on='node_id')
        edges_cluster = edges_cluster.drop('node_id',axis=1).rename(columns={cluster_type : 'start_comp'})
        
        edges_cluster = pd.merge(edges_cluster, nodes_gdf[['node_id',cluster_type]],left_on='ds_node_id',right_on='node_id')
        edges_cluster = edges_cluster.drop('node_id',axis=1).rename(columns={cluster_type : 'end_comp'})
                
        ind = edges_cluster.end_comp == edges_cluster.start_comp
        
        edges_cluster.loc[ind, 'compart_id'] = edges_cluster.loc[ind, 'start_comp']
        edges_between = edges_cluster.loc[~ind]
        
        """Aggregate nodes to compartment
        """
        def f(x):
            d = {}
            d['node_stor'] = x['storage'].sum()
            d['initial_storage'] = x['initial_storage'].sum()
            d['surf_area'] = x['surf_area'].sum()
            d['run_coef'] = (x['surf_area'] * x['run_coef'] / d['surf_area']).sum()
            d['geometry'] = Point((x.geometry.x).mean(),(x.geometry.y).mean())
            d['chamber_ar'] = x['chamber_ar'].sum()
            
            return pd.Series(d)
        compartment_nodes = nodes_gdf.groupby(cluster_type).apply(f).reset_index().rename(columns={cluster_type : 'compart_id'})

        #Add in pipe storage
        pipe_stor = edges_cluster.groupby('compart_id').sum().storage.reset_index()
        compartment_nodes = pd.merge(compartment_nodes, pipe_stor.rename(columns={'storage':'pipe_stor'}), on = 'compart_id', how="outer")
        compartment_nodes.pipe_stor = compartment_nodes.pipe_stor.fillna(0)
        compartment_nodes['storage'] = compartment_nodes[['node_stor','pipe_stor']].sum(axis=1)
        del pipe_stor
        
        """Calculate time area graph for each compartment
        """
        compartment_nodes['timearea'] = None
        compartment_nodes['timearea_pipe'] = None
        compartment_nodes['timearea_surface'] = None
        compartment_nodes['pipe_time'] = 0
        def update_timearea(timeareas, cumul_pipe_travel, outfall, nodes, edges):
            #Recursive function to trace upstream and create timearea graph
            
            node = outfall.us_node_id
            row = outfall
            if node not in timeareas.keys():
                
                cumul_pipe_travel += row.travel_time # s
                
                timeareas[node] = {}
                timeareas[node]['pipe_time'] = cumul_pipe_travel
                timeareas[node]['surface_time'] = SURFACE_TRAVEL_TIME
                timeareas[node]['surf_area'] = nodes.loc[node,['run_coef','surf_area']].prod()
                outfalls = edges.loc[edges.ds_node_id == node]
                for idx, edge in outfalls.iterrows():
                    update_timearea(timeareas, cumul_pipe_travel, edge, nodes, edges)
                
                
        for compart_id in compartment_nodes.compart_id:
            in_elev = None
            out_elev = None
            nodes_in_compartment = nodes_gdf.loc[nodes_gdf[cluster_type] == compart_id].set_index('node_id')
            edges_in_compartment = edges_cluster.loc[edges_cluster.start_comp == compart_id]
            edges_in_compartment_ = edges_cluster.loc[edges_cluster.end_comp == compart_id]
            if edges_in_compartment.size > 0:
                outfalls = edges_between.loc[edges_between.info_id.isin(edges_in_compartment.info_id)].copy()
                if out_elev is None:
                    out_elev = nodes_in_compartment.loc[outfalls.us_node_id, 'chamber_fl'].mean()
                    compartment_nodes.loc[compart_id,'out_elev'] = out_elev
                outfalls_ = edges_between.loc[edges_between.info_id.isin(edges_in_compartment_.info_id)].copy()
                if in_elev is None:
                    in_elev = nodes_in_compartment.loc[outfalls_.ds_node_id,'chamber_fl'].mean()
                    compartment_nodes.loc[compart_id,'in_elev'] = in_elev
                outfalls.travel_time = 0 #Ignore edges between travel time (since they will be modelled live)
                timeareas = {}
                for idx, edge in outfalls.iterrows():
                    cumul_pipe_travel = 0
                    update_timearea(timeareas, cumul_pipe_travel, edge, nodes_in_compartment, edges_in_compartment)
                
                timeareas = pd.DataFrame(timeareas).T
                
                times = linspace(0,timeareas.pipe_time.max() + timeareas.surface_time.max(),TIMEAREA_RESOLUTION)
                areas = zeros(times.shape)
                
   
                """Pipe flow timearea diagram
                """
                #Superimpose timearea diagram
                for idx, row in timeareas.iterrows():
                    ind_start = where(times >= row.pipe_time)[0][0]
                    ind_end = where(times >= (row.pipe_time + row.surface_time))[0][0]
                    areas[ind_start : ind_end] += linspace(0,row.surf_area, ind_end - ind_start)
                    areas[ind_end:] += row.surf_area

                #Calculate timearea diagram for entry to pipe network
                surface_areas = zeros(times.shape)
                
                for time, area in timeareas.sort_values(by='surface_time').set_index('surface_time').surf_area.iteritems():
                    ind_start = 0
                    ind_end = where(times >= time)[0][0]
                    surface_areas[ind_start : ind_end] += linspace(0,area, ind_end - ind_start)
                    surface_areas[ind_end:] += area
                    
                def timearea_format(t, a):
                    a_ = diff(a)
                    ta = pd.DataFrame(a_, index=pd.to_datetime(t[1:], unit='s'), columns = ['surf_area'])
                    ta = ta.resample(DT_SIM).sum()
                    if ta.surf_area.sum() > 0:
                        ta['normalised'] = ta.surf_area.div(ta.surf_area.sum())
                    else:
                        ta['normalised'] = 0
                    ta.index = range(ta.shape[0])
                    ta = ta.loc[ta.normalised > 0]
                    return ta
                
                timearea_surface_runoff = timearea_format(times, surface_areas)
                timearea_outflow = timearea_format(times, areas)
                
                #Warning - this is approximate!
                timearea_pipe_flow = timeareas.copy()
                timearea_pipe_flow.index = pd.to_datetime(timearea_pipe_flow.pipe_time,unit='s')
                timearea_pipe_flow.index += DT_SIM/2
                timearea_pipe_flow = timearea_pipe_flow.resample(DT_SIM).sum()
                timearea_pipe_flow = timearea_pipe_flow.loc[timearea_pipe_flow.surf_area > 0]
                timearea_pipe_flow['normalised'] = timearea_pipe_flow.surf_area.div(timearea_pipe_flow.surf_area.sum())
                timearea_pipe_flow.index = range(timearea_pipe_flow.shape[0])
                
                compart_index = compartment_nodes.loc[compartment_nodes.compart_id == compart_id].index[0]
                compartment_nodes.at[compart_index, "timearea"] = timearea_outflow.normalised.to_dict() #Use a dict because can be stored in geojson
                compartment_nodes.at[compart_index, "timearea_pipe"] = timearea_pipe_flow.normalised.to_dict() #Use a dict because can be stored in geojson
                compartment_nodes.at[compart_index, "timearea_surface"] = timearea_surface_runoff.normalised.to_dict() #Use a dict because can be stored in geojson
                entry_nodes = edges_cluster.loc[(edges_cluster.end_comp == compart_id) & (edges_cluster.start_comp != compart_id), 'ds_node_id'].values
                entry_nodes = timeareas.index.intersection(entry_nodes) #Only count time from nodes that are physically connected to outfalls
                if entry_nodes.size > 0:
                    compartment_nodes.at[compart_index, "pipe_time"] = timeareas.loc[entry_nodes, 'pipe_time'].mean() * S_TO_DT
            else:
                compartment_nodes.loc[compart_id, "in_elev"] = nodes_in_compartment.chamber_fl.mean()
                compartment_nodes.loc[compart_id, "out_elev"] = nodes_in_compartment.chamber_fl.mean()
        """Define pipe network
        """
        starts = compartment_nodes.loc[edges_between.start_comp].geometry
        ends = compartment_nodes.loc[edges_between.end_comp].geometry
        geom = [LineString([x,y]) for x,y in zip(starts,ends)]
        
        edges_between = gpd.GeoDataFrame(data = edges_between.values, 
                                         columns = edges_between.columns, 
                                         geometry = geom, 
                                         index = edges_between.index) # Don't judge me... I get a slice warning if I do this any other way
    
        edges_between['number_of_timesteps'] = round(edges_gdf.loc[edges_between.index,'travel_time'] * S_TO_DT)
        
        """Calculate upstreamness
        """
        #Assign upstreamness
        outfall_ind = nodes_gdf.loc[nodes_gdf.node_type == 'Outfall',cluster_type].values
        upstreamness = {x : 0 for x in compartment_nodes.loc[outfall_ind].index.tolist()}
        
        def assign_upstream(links, upstreamness):
            upstreamness_ = upstreamness.copy()
            in_nodes = links.loc[links['end_comp'].isin(upstreamness.keys()),'start_comp']
            in_nodes = list(set(in_nodes).difference(upstreamness.keys()))
            ind = max(list(upstreamness_.values())) + 1
            for node in in_nodes:
                upstreamness[node] = ind
            if upstreamness == upstreamness_:
                return upstreamness
            else:
                upstreamness = assign_upstream(links,upstreamness)
                return upstreamness
        
        upstreamness = assign_upstream(edges_between,upstreamness)
        compartment_nodes['upstr_ind'] = pd.Series(upstreamness)
        del upstreamness, outfall_ind
        
        
        fid = __file__
        fn = os.path.basename(fid)
        copyfile(fid, os.path.join(cluster_folder, fn))
        
        """Write Compartments
        """
        gpd.GeoDataFrame(compartment_nodes, crs=CRS).to_file(driver = driver, filename= os.path.join(cluster_folder,"compartment_nodes" + extension))
        gpd.GeoDataFrame(edges_between, crs=CRS).to_file(driver = driver, filename= os.path.join(cluster_folder,"compartment_edges" + extension))
                
        #Front bit of this could definitely be done out of the loop
        """Aggregate detailed results to compartments
        """
    
        for rain in sim_rains:
            rain_folder = os.path.join(cluster_folder, rain)
            if not os.path.exists(rain_folder):
                os.mkdir(rain_folder)
            
            #I've removed 'floodloss' for now because it's constantly 0
            info_node = {'floodvol' : os.path.join(raw_root,info_dir, rain,"floodvolume.csv"),
                              'depth' : os.path.join(raw_root,info_dir, rain, "depnod.csv"),
                              'volume' : os.path.join(raw_root,info_dir, rain, "volume.csv")}
            info_flow = os.path.join(raw_root,info_dir,rain,"ds_flow.csv")
            
            def load_info(fid, type_):
                temp_df = pd.read_csv(fid)
                try:
                    temp_df.index = pd.to_datetime(temp_df.Time, format = "%d/%m/%Y %H:%M:%S")
                except:
                    temp_df.index = pd.to_datetime(temp_df.Time.replace({'00/00/0000' : sim_start,
                                                                         '01/00/0000' : sim_start[0:8] + '02'},
                                                                        regex=True))
                temp_df.index.name = 'time'
                temp_df = temp_df.drop(['Time', 'Seconds'], axis=1)
                
                temp_df = temp_df.resample(DT_RAW).mean()
                if type_ == 'arc':
                    temp_df = temp_df.mul(DT_RAW / DT_ARC) # convert from mean over DT to total over DT
                    if DT_SIM < DT_RAW:
                        temp_df = temp_df.resample(DT_SIM,label='right').interpolate(limit = int(DT_RAW/DT_SIM))*(DT_SIM/DT_RAW)
                    else:
                        temp_df = temp_df.resample(DT_SIM,label='right').sum(min_count = 1)
                elif type_ == 'node_id':
                    if DT_SIM < DT_RAW:
                        temp_df = temp_df.resample(DT_SIM,label='right').interpolate(limit = int(DT_RAW/DT_SIM))
                    else:
                        temp_df = temp_df.resample(DT_SIM,label='right').mean()
                
                # temp_df = temp_df.astype(float16) # unstacking the larger datasets is pretty heavyweight

                return temp_df
            
            #Load and convert flow data
            flow_df = load_info(info_flow,'arc')
            flow_df = flow_df[flow_df.columns.intersection(edges_gdf.info_id)]
            
            compart_flows = flow_df[flow_df.columns.intersection(edges_between.info_id)].unstack().reset_index().rename(columns={'level_0' : 'arc', 0 : 'val', 'Seconds' : 'time'})
        
            #Write flows between compartments
            compart_flows = pd.merge(edges_between[['start_comp','end_comp','info_id']], compart_flows, left_on ="info_id", right_on="arc")
            
            #Remove between compartment flows and add to node volumes
            internal_flows = flow_df[flow_df.columns[~flow_df.columns.isin(compart_flows.info_id.unique())]]
            internal_nodes = edges_gdf.loc[edges_gdf.info_id.isin(internal_flows.columns), 'us_node_id'].unique()
            node_flows = pd.DataFrame(index = internal_flows.index, 
                                      columns = internal_nodes,
                                      data = zeros((internal_flows.shape[0], len(internal_nodes))))
            
            for arc in internal_flows.columns.intersection(flow_df.columns):
                cap = edges_gdf.loc[edges_gdf.info_id == arc,'capacity'].values[0]
                full_storage = edges_gdf.loc[edges_gdf.info_id == arc,['length','cross_sec']].prod(axis=1).values[0]
                if cap == 0:
                    pct_full = pd.Series(index = internal_flows.index,
                                         data = zeros(internal_flows.shape[0]))
                    pct_full.loc[internal_flows[arc] > cap] = 1
                else:
                    pct_full = internal_flows[arc] / cap
                    pct_full.loc[pct_full > 1] = 1
                pipe_s = (pct_full * full_storage).abs()
                us_node = edges_gdf.loc[edges_gdf.info_id == arc,'us_node_id'].values[0]
                node_flows[us_node] += pipe_s

            compartments = nodes_gdf[cluster_type].unique()
            compartment_nodes = []
            for key in info_node.keys():
                node_df = load_info(info_node[key],'node_id')
                if key == 'volume':
                    node_df = node_df.add(node_flows, fill_value = 0)
                elif (key == 'floodvol') | (key == 'depth'):
                    #if floodvol < 0 then that means avaiable capacity
                    #no idea what negative depth means... best avoided
                    node_df = node_df.clip(lower = 0) 
                    
                compart_node = pd.DataFrame(columns = compartments, 
                                         index = node_flows.index, 
                                         data = zeros((node_flows.shape[0], len(compartments))))
                compart_node[:] = NaN
                for cid in compartments:
                    nodes = nodes_gdf.loc[nodes_gdf[cluster_type] == cid, 'node_id'].values
                    nodes = set(nodes).intersection(node_df.columns)
                    if key in ['volume','floodvol']:
                        compart_node[cid] = node_df[nodes].sum(axis=1, min_count=1)
                    elif key == 'depth':
                        compart_node[cid] = node_df[nodes].mean(axis=1)
                    
                compart_node = compart_node.unstack().reset_index().rename(columns={'level_0' : cluster_type, 0 : key})
                compartment_nodes.append(compart_node)
            compartment_nodes_ = pd.concat(compartment_nodes, axis=1)
            compartment_nodes_ = compartment_nodes_.loc[:,~compartment_nodes_.columns.duplicated()]
       
            """Write aggregated results
            """
            comparison_arcs = pd.read_csv(os.path.join(cluster_root, "comparison_arcs.csv"))
            comparison_nodes = pd.read_csv(os.path.join(cluster_root, "comparison_nodes.csv"))
            
            compart_flows = compart_flows.loc[compart_flows.arc.isin(comparison_arcs.info_id)]
            compart_flows.to_parquet(os.path.join(rain_folder,"highfid_flows.gzip"), index=True, compression='GZIP')
            
            node_lookup = nodes_gdf[[cluster_type,'node_id']].drop_duplicates()
            node_lookup = node_lookup.loc[node_lookup.node_id.isin(comparison_nodes.node_id)]
            
            compartment_nodes_ = compartment_nodes_.loc[compartment_nodes_[cluster_type].isin(node_lookup[cluster_type])]
            compartment_nodes_.to_parquet(os.path.join(rain_folder,"highfid_nodes.gzip"), index=True, compression='GZIP')
       
 
if __name__ == "__main__":
    
    if len(sys.argv) > 1:
        runid = str(sys.argv[1])
        jobid = int(sys.argv[2])
        njobs = int(sys.argv[3])
        submi = int(sys.argv[4])
        submimax = int(sys.argv[5])
    else:
        runid = '2021-03-18'
        jobid = 0
        njobs = 1
        submi = 0
        submimax = 1
    
    root = os.path.join("/rds","general","user","bdobson", "home", "cwsd_sewer","data")
    catchment = "cranbrook"
    sim_root = os.path.join(root,catchment,"results",runid)
    dts = [5,10,20,50,100,300,1440,2,1,0.5]
    # dts = [1000,1440]
    clusters = sorted(glob(os.path.join(sim_root,"*","*.geojson")))
 
    jobs = []
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
                  jobs.append({'dt' : dt,
                               'cluster' : col,
                               'root' : os.path.split(cluster)[0],
  			     'jobid' : str(jobid)})
    
    jobs = jobs[int(submi/submimax * len(jobs)):int((submi+1)/submimax * len(jobs))]    
    for (i, job) in enumerate(jobs):
        if (i % njobs) == jobid:
            print('starting')
            print(job)
            
            aggregate(root, job['root'], job['cluster'], catchment, dt_sim = job['dt'])
            print('ending')
            print(job)
    
    
