# -*- coding: utf-8 -*-
"""
Created on Fri Jan 22 08:08:18 2021

@author: bdobson
"""
import constants
from random import shuffle
from numpy import isnan
from tqdm import tqdm
#Contains all code required to simulate a partitioned CWSD sewer network
class model:
    """Class that contains cwsd nodes and arcs with functions to add them from generic databases
    """
    def __init__(self, physical = False):
        self.model_arcs = {}
        self.model_nodes = {}
        self.model_nodes_type = {}
        for node in Node.__subclasses__():
            self.model_nodes_type[node.__name__] = {}
        
        self.physical = physical
            
        
    def add_nodes(self, nodes_dict):
        self.upstr_dict = {}
        for name, data in nodes_dict.items():
                
            
            name = str(name)
            if data['surf_area'] > 0:
                node_name = name + '-land'
                self.model_nodes_type['Land'][node_name] = Land(**data)
                self.model_nodes[node_name] = self.model_nodes_type['Land'][node_name]
                self.model_nodes[node_name].timearea_surface = {int(k):v for k,v in self.model_nodes[node_name].timearea_surface.items()}
                
                ta_sum = sum(self.model_nodes[node_name].timearea_surface.values())
                for time, entry in self.model_nodes[node_name].timearea_surface.items():
                    self.model_nodes[node_name].timearea_surface[time] = entry / ta_sum
                
                self.model_nodes[node_name].name = node_name
                self.model_nodes[node_name].data_input = 'rain_holder'
                
            if data['compartment_type'] == 'basic':
                node_name = name + '-sewer'

                if not data['timearea_pipe']:
                    data['timearea_pipe'] = {0:1}
                self.model_nodes_type['Sewerage'][node_name] = Sewerage(**data)
                self.model_nodes[node_name] = self.model_nodes_type['Sewerage'][node_name]
                self.model_nodes[node_name].name = node_name
                self.model_nodes[node_name].timearea_pipe = {int(k):v for k,v in self.model_nodes[node_name].timearea_pipe.items()}
                
                ta_sum = sum(self.model_nodes[node_name].timearea_pipe.values())
                for time, entry in self.model_nodes[node_name].timearea_pipe.items():
                    self.model_nodes[node_name].timearea_pipe[time] = entry / ta_sum
                
                upstr_ind = self.model_nodes[node_name].upstr_ind
                if upstr_ind in self.upstr_dict.keys():
                    self.upstr_dict[upstr_ind].append(self.model_nodes[node_name])
                else:
                    self.upstr_dict[upstr_ind] = [self.model_nodes[node_name]]
                self.upstr_dict_keys = list(self.upstr_dict.keys())
                self.upstr_dict_keys.sort(reverse = True) # True is downstream to upstream, False is upstream to downstream
                
            elif data['compartment_type'] == 'outfall':
                node_name = name + '-outfall'
                self.model_nodes_type['Junction'][node_name] = Junction(**data)
                self.model_nodes[node_name] = self.model_nodes_type['Junction'][node_name]
                self.model_nodes[node_name].name = node_name
                
    def add_arcs(self, arcs_dict):
        outfalls = [int(x.replace('-outfall','')) for x in self.model_nodes_type['Junction'].keys()]
        #Define network
        for name, arc in arcs_dict.items():
            s_Arc = Arc
            nAR = None
            if arc['end_comp'] in outfalls:
                name_append = '-outfall'
                
            else:
                name_append = '-sewer'
                if self.physical:
                    s_Arc = Arc_physical
                    A = arc['cross_sec']
                    r = (A / constants.PI) ** 0.5
                    P = r * constants.PI * 2
                    if P <= 0:
                        if arc['gradient'] > 0:
                            nAR = (arc['capacity'] / constants.M3_S_TO_M3_DT) / (arc['gradient'] ** 0.5) # If looks dodgy, then back calculate
                        else:
                            nAR = (arc['capacity'] / constants.M3_S_TO_M3_DT) / (0.001 ** 0.5) # generic slope
                    else:
                        nAR = arc['mannings'] * A * ((A/P) ** (2/3))
                    if arc['edge_type'] == 'weir':
                        s_Arc = Arc_weir
                    elif arc['edge_type'] == "river_reach_link":
                        s_Arc = Arc
                    
            self.model_arcs[name] = s_Arc(name=name,
                                            inPort=self.model_nodes[str(arc['start_comp']) + '-sewer'],
                                            outPort=self.model_nodes[str(arc['end_comp']) + name_append],
                                            capacity=arc['capacity'],
                                            preference=1,
                                            cross_section = arc['cross_sec'],
                                            length = arc['length'],
                                            number_of_timesteps = int(arc['number_of_timesteps']),
                                            gradient = arc['gradient'],
                                            mannings = arc['mannings'],
                                            mannings_coef_nAR = nAR,
                                            type_ = arc['edge_type'],
                                            # volumetric_capacity = arc['storage'] # THERE ARE ARCS WITH 0 HERE... OBVIOUSLY AN ISSUE
                                            )
            if (arc['edge_type'] == 'weir') & self.physical:
                self.model_arcs[name].weir_param_mul = float(arc['weir_param_mul'])
                self.model_arcs[name].weir_param_pow = float(arc['weir_param_pow'])
                self.model_arcs[name].weir_height = float(arc['weir_height'])
                s_Arc = Arc_physical
            elif arc['edge_type'] == "river_reach_link":
                s_Arc = Arc_physical
                
        #Land to sewer and sewer to land
        for key, land_node in self.model_nodes_type['Land'].items():
            compartment = land_node.name.replace('-land','')
            name = compartment + '-runoff'
            self.model_arcs[name] = Arc(name=name,
                                        inPort = land_node,
                                        outPort = self.model_nodes[compartment + '-sewer'],
                                        capacity = constants.UNBOUNDED_CAPACITY,
                                        number_of_timesteps = 0,
                                        preference = 1,
                                        type_ = 'runoff')
            name = compartment + '-spill'
            self.model_arcs[name] = Arc(name=name,
                                        inPort = self.model_nodes[compartment + '-sewer'],
                                        outPort = land_node,
                                        capacity = constants.UNBOUNDED_CAPACITY,
                                        number_of_timesteps = 0,
                                        preference = 1,
                                        type_ = 'spill')
            
    def add_inputs(self, inputs_dict):
        dates = []
        for node in self.model_nodes.values():
            if hasattr(node,'data_input'):
                node.data_input_dict = inputs_dict[node.data_input]
                dates.append(node.data_input_dict.keys())
                        
        #Inputs must be defined for all dates
        self.dates = list(set.intersection(*[set(x) for x in dates]))
        
        #Relies on dates being in ISO 8601 to sort correctly
        self.dates.sort()

    def run(self):
        flows = []
        storages = []
        flood_vol = []
        depth = []
        for date in tqdm(self.dates):
            #Update time
            for node in self.model_nodes.values():
                node.date = date

            #Create runoff
            for name, node in self.model_nodes_type['Land'].items():
                node.create_runoff()
                
            
            #Discharge sewers in upstreamness order
            for key in self.upstr_dict_keys:
                nl = list(self.upstr_dict[key])
                shuffle(nl)

                for node in nl:
                    node.make_discharge()
                    
            
            for name, arc in self.model_arcs.items():
                flows.append({'date' : date,
                              'arc' : name,
                              'flow_in' : arc.volume_in,
                              'flow_out' : arc.volume_out})
            #Mass balance
            running_in = 0
            running_out = 0
            running_ds = 0
            for name, node in self.model_nodes.items():
                if 'outfall' not in name:
                    totin = 0 
                    totout = 0
                    ds = 0
                    for arc in node.inArcs.values():
                        totin += arc.volume_out
                    for arc in node.outArcs.values():
                        totout += arc.volume_in
                    
                    if node.type_ == "Land":
                        ds = (node.storage['volume'] - node.storage_['volume'])
                        flood_vol.append({'date' : date, 
                                          'node' : name,
                                          'val' : node.storage['volume']})
                        totin += (node.rain['volume'] - node.greenspace_rain)
                    elif node.type_ == "Sewerage":
                        ds = (node.storage['volume'] - node.storage_['volume'])
                        depth.append({'date' : date, 
                                          'node' : name,
                                          'val' : node.get_head('av')})
                        storages.append({'date' : date,
                                     'node' : node.name,
                                     'val' : node.storage['volume']})
                    if (totin - totout - ds) > constants.FLOAT_ACCURACY:
                        print("Mass balance error at {} on {}".format(name, date))
                    
                    running_in += totin
                    running_out += totout
                    running_ds += ds

            if (running_in - running_out - running_ds) > constants.FLOAT_ACCURACY:
                print('system mass balance error at date : {}'.format(date))

            #end timesteps
            for node in self.model_nodes.values():
                node.end_timestep()
            
            for arc in self.model_arcs.values():
                arc.end_timestep()

        return {'flows' : flows,
                'storages' : storages,
                'flood_vol' : flood_vol,
                'depths' : depth}
class Node:
    def __init__(self,**kwargs):
        self.inArcs = {}
        self.outArcs = {}
        self.name = None
        self.date = None
        self.flag = 0
        self.query_value = 0
        self.month = 0
        self.year = 0
        self.losses = 0
        self.__dict__.update(kwargs)

    def read_input(self):
        self.query_value = self.data_input_dict[self.date]
        return self.query_value
    
    def get_connected(self,of_type = None):
        priorities = {}
        total_avail = 0
        total_priority = 0
            
        for name, arc in self.inArcs.items():
            if (arc.inPort.type == of_type) | (of_type is None):
                avail = arc.send_pull_check()
                if arc.preference < constants.FLOAT_ACCURACY:
                    avail = 0
                priorities[name] = avail*arc.preference
                total_avail += avail
                total_priority += avail*arc.preference
        
        return {'total_avail' : total_avail, 'total_priority' : total_priority,'priorities' : priorities}

    def end_timestep(self):
        pass

    def raw_concentration(self):
        return {'volume' : 0,
                'misc_pol' : 0}

    def update_concentration(self,volume):
        concentration = self.raw_concentration()
        concentration['volume'] = volume
        return concentration
    
    def copy_concentration(self,c):
        return dict([('volume',c['volume'])] + [(key,c[key]) for key in constants.POLLUTANTS])
    
    def empty_concentration(self):
        return dict([('volume',0)] + [(key,0) for key in constants.POLLUTANTS])
    
    def blend_concentrations(self, c1, c2):
        c = self.empty_concentration()
        
        c['volume'] = c1['volume'] + c2['volume']
        if c['volume'] > constants.FLOAT_ACCURACY:
            for pollutant in constants.POLLUTANTS:
                c[pollutant] = (c1[pollutant]*c1['volume'] + c2[pollutant] * c2['volume'])/c['volume']
            
        return c
    
class Sewerage(Node):
    def __init__(self,**kwargs):
        self.leakage = 0
        self.discharge_preference_type = 'proportional' #Can be 'absolute' or 'proportional'
        self.type_ = 'Sewerage'
        self.storage = 0
        self.pipe_stor = 0
        self.initial_storage = 0
        super().__init__(**kwargs)
        self.storage_ = self.update_concentration(self.initial_storage)
        self.storage = self.copy_concentration(self.storage_)
        self.pipe_queue_length = max([max(map(int,self.timearea_pipe.keys())) + 1, 
                                      self.pipe_time + 1,
                                      2])
        self.storage_queue = [ 0 for v in range(self.pipe_queue_length)]
        self.storage_queue[0] = self.storage['volume']
        
    def end_timestep(self):
        self.storage_ = self.copy_concentration(self.storage)
        self.storage_queue[1] += self.storage_queue[0]
        self.storage_queue = self.storage_queue[1:] + [0]
    
    def set_push_request(self, sender, request):

        request_ = self.copy_concentration(request)
        
        
        self.storage = self.blend_concentrations(self.storage, request_)

        
        if sender.inPort.type_ == "Land":

            for time, normalised in self.timearea_pipe.items():
                self.storage_queue[time] += request_['volume'] * normalised
        elif sender.inPort.type_ == "Sewerage":
            self.storage_queue[max(self.pipe_time,0)] += request_['volume']                
        
        return self.empty_concentration()
    
    def make_discharge(self):
        if (sum(self.storage_queue) - self.storage['volume']) > constants.FLOAT_ACCURACY:
            print('storage_queue doesnt match storage for ' + self.name)
        
        discharge = self.copy_concentration(self.storage)
        discharge['volume'] = self.storage_queue[0]
        
        rejected = self.copy_concentration(discharge)
        rejected['volume'] = 0
        
        preference = {}
        for name, arc in self.outArcs.items():
            if 'spill' not in name:
                preference[name] = arc.getExcess()
        
        push_capacity = sum(preference.values())
        push_actual = min(discharge['volume'], push_capacity)
        if push_actual < 0:
            print(" ".join(["Attempted negative push at",self.name,self.date]))
            push_actual = 0
        
        if push_actual > 0:
            for name, arc in self.outArcs.items():
                if 'spill' not in name:
                    discharge_ = self.copy_concentration(discharge)
                    discharge_['volume'] = push_actual * preference[name] / push_capacity
                    reply = arc.send_push_request(discharge_)
                    rejected['volume'] += reply['volume']
            
            discharge['volume'] -= (push_actual - rejected['volume'])
            self.storage['volume'] -= (push_actual - rejected['volume'])
                
        
        #Spill leftover
        discharge_ = self.copy_concentration(discharge)
        discharge_['volume'] = max(self.storage['volume'] - self.capacity,0)
        discharge_['volume'] = min(discharge['volume'], discharge_['volume'])
        spill_vol = discharge_['volume']
        if discharge_['volume'] > 0:
            for name, arc in self.outArcs.items():
                if 'spill' in name:
                    #This assumes one spill only
                    reply = arc.send_push_request(discharge_)
                    discharge['volume'] -= (spill_vol - reply['volume'])
                    self.storage['volume'] -= (spill_vol - reply['volume'])
        self.storage_queue[0] = discharge['volume']
          
    def get_head(self,typeof):
    
        
        if typeof == 'in':
            elev = self.in_elev
        elif typeof == 'out':
            elev = self.out_elev
        elif typeof == 'av':
            if isnan(self.in_elev):
                if isnan(self.out_elev):
                    print('no elevation')
                else:
                    elev = self.out_elev
            elif isnan(self.out_elev):
                elev = self.in_elev
            else:
                elev = (self.in_elev + self.out_elev)/2
        else:
            print("Warning: specify typeof gethead()")
        

        if self.storage['volume'] > self.pipe_stor:
            elev += (self.storage['volume'] - self.pipe_stor) / self.chamber_ar

            return elev
        else:
            return elev
            
class Land(Node):
    def __init__(self,**kwargs):
        self.type_ = 'Land'
        self.surf_area = 0
        self.run_coef = 0
        self.timearea = [{0:0}]
        super().__init__(**kwargs)
        self.storage = self.empty_concentration()
        self.storage_ = self.empty_concentration()
        self.greenspace_rain = 0
        self.capacity = constants.UNBOUNDED_CAPACITY
    def get_push_available(self):
        avail = self.capacity - self.storage['volume']
        return avail
    
    def create_runoff(self):
          
        #Clear storage
        for arc in self.outArcs.values():
            storage = self.copy_concentration(self.storage)
            storage['volume'] *= arc.preference
            
            self.storage = arc.send_push_request(storage, 0)

        #Read rainfall
        conversion = self.surf_area * constants.MM_M2_TO_SIM_VOLUME
        self.rain = self.update_concentration(conversion * self.read_input())
        
        #Send rain to sewer
        runoff_to_sewer = self.rain['volume'] * self.run_coef
        runoff_to_sewer = self.update_concentration(runoff_to_sewer)    
        
        for arc in self.outArcs.values():
            runoff_to_sewer_ = self.copy_concentration(runoff_to_sewer)
            runoff_to_sewer_['volume'] *= arc.preference # Send along arcs according to preference
            
            for time, normalised in self.timearea_surface.items():
                temp = self.copy_concentration(runoff_to_sewer_)
                temp['volume'] = runoff_to_sewer_['volume'] * normalised
                reply = arc.send_push_request(temp, time) # This will cause sewer storage to be underestimated, because in-pipe travel time is included in the timearea time
                self.storage = self.blend_concentrations(reply, self.storage)
                        
        #Track greenspace rain
        self.greenspace_rain = self.rain['volume'] * (1 - self.run_coef)
    
    def set_push_request(self, sender, request):
        self.storage = self.blend_concentrations(request, self.storage)
        return self.empty_concentration()
    
    def end_timestep(self):
        self.storage_ = self.copy_concentration(self.storage)        
        

class Junction(Node):
    def __init__(self,**kwargs):
        self.type_ = 'Junction'
        super().__init__(**kwargs)
        
        self.storage = self.empty_concentration()
        
        self.capacity = constants.UNBOUNDED_CAPACITY
        
    def get_push_available(self):
        return constants.UNBOUNDED_CAPACITY
    
    def set_push_request(self, sender, request):
        #Push requests take concentrations as an argument
        if 'outfall' in self.name:
            spill = self.empty_concentration()
        
        return spill
    
class Arc:
    def __init__(self,**kwargs):
        self.number_of_timesteps = 0
        self.queued_pushes = []
        self.flow_out = 0
        self.flow_in = 0
        self.volume_in = 0
        self.volume_out = 0
        
        self.length = 0
        self.cross_section = 0
        self.__dict__.update(kwargs)
        self.m_capacity = self.capacity
        self.inPort.outArcs[self.name] = self
        self.outPort.inArcs[self.name] = self
            
    def getExcess(self):
        excess = self.capacity - self.flow_in
        outPort_expected_push = 0
        for name, arc in self.outPort.outArcs.items():
            if 'spill' not in name:
                outPort_expected_push += arc.capacity
        
        excess = min(excess, max(self.outPort.capacity + outPort_expected_push - self.outPort.storage['volume'],0))

        return excess
    
    def update_queue(self, push):
        push['average_flow'] = push['request']['volume'] / (push['timesteps'] + 1)
        self.flow_in += push['average_flow']
        self.volume_in += push['request']['volume'] 
        self.queued_pushes.append(push)
        
    def copy_concentration(self,c):
        return dict([('volume',c['volume'])] + [(key,c[key]) for key in constants.POLLUTANTS])

    
    def empty_concentration(self):
        return dict([('volume',0)] + [(key,0) for key in constants.POLLUTANTS])
       
    
    def send_push_request(self, request, timesteps = 0):

        excess_in = self.getExcess()
        
        not_pushed = self.copy_concentration(request)
        not_pushed['volume'] = max(0, request['volume'] - excess_in)
        
        request['volume'] -= not_pushed['volume']
        
        self.update_queue({'timesteps' : timesteps + self.number_of_timesteps, 'request' : request})
        done_pushes = []
        for push in self.queued_pushes:
            if push['request']['volume'] < constants.FLOAT_ACCURACY:
                done_pushes.append(push)
            elif push['timesteps'] == 0:
                excess_out = self.m_capacity - self.flow_out
                original_request = self.copy_concentration(push['request'])
                push_volume = min(original_request['volume'], excess_out)
                backup_volume = original_request['volume'] - push_volume
                push['request']['volume'] -= backup_volume
                reply = self.outPort.set_push_request(self, push['request'])
                
                not_pushed['volume'] += (backup_volume + reply['volume'])
                
                push['request']['volume'] = 0
                total_removed = original_request['volume'] - (reply['volume'] + backup_volume)
                self.flow_out += (push['average_flow'] * total_removed / original_request['volume'])
                self.volume_out += total_removed

        for push in done_pushes:
            self.queued_pushes.remove(push)
        return not_pushed
    
    def end_timestep(self):

        self.flow_in = 0
        self.flow_out = 0
        self.volume_in = 0
        self.volume_out = 0
        for push in self.queued_pushes:
            push['timesteps'] = max(push['timesteps'] - 1, 0)

class Arc_physical(Arc):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        
        
    def getExcess(self):
        inhead = self.inPort.get_head('out')
        outhead = self.outPort.get_head('in')
        slope = (inhead - outhead)/self.length

        if slope <= 0:
            slope = self.gradient
        
        
        self.m_capacity = (self.mannings_coef_nAR * (slope ** 0.5)) * constants.M3_S_TO_M3_DT
        
        excess = self.m_capacity - self.flow_in
        
        outPort_expected_push = 0
        for name, arc in self.outPort.outArcs.items():
            if 'spill' not in name:
                outPort_expected_push += arc.capacity
        
        excess = min(excess, max(self.outPort.capacity + outPort_expected_push - self.outPort.storage['volume'],0))
        
        return excess

class Arc_weir(Arc):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        
        
    def getExcess(self):
        
        
        
        elev = max(self.inPort.get_head('out') - self.weir_height, 0)

        self.m_capacity = self.weir_param_mul * (elev ** self.weir_param_pow) * constants.M3_S_TO_M3_DT
        
        excess = self.m_capacity - self.flow_in
        outPort_expected_push = 0
        for name, arc in self.outPort.outArcs.items():
            if 'spill' not in name:
                outPort_expected_push += arc.capacity
        
        excess = min(excess, max(self.outPort.capacity + outPort_expected_push - self.outPort.storage['volume'],0))
        
        return excess
