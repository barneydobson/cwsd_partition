# -*- coding: utf-8 -*-
"""
Created on Thu Mar 18 08:53:07 2021

@author: bdobson
"""
import numpy as np
from sknetwork.clustering import Louvain
from sklearn import cluster
import networkx as nx
import warnings
#Implementation of various partitioning algorithms

def wrapper(adj_mat,
             n,
             coords,
             method):
    
    params = {'n' : n,
              'adj_mat' : adj_mat,
              'coords' : coords}
    
    #Set functions and params
    if method == 'louv':
        cf = louvain
    elif method == 'sc':
        cf = spectral        
    elif method == 'km':
        cf = geo_cluster
        params['method'] = 'KMeans'
    elif method == 'ac':
        cf = geo_cluster
        params['method'] = 'AgglomerativeClustering'
    elif method == 'L_L_sc':
        cf = L_L_sc
    if method == 'L_scn':
        cf = L_scn
    elif method == 'km_sc':
        cf = geo_cluster_sc
        params['method'] = 'KMeans'
    elif method == 'ac_sc':
        cf = geo_cluster_sc
        params['method'] = 'AgglomerativeClustering'
    elif method == 'fluid':
        cf = fluid

    warnings.filterwarnings(
                    "ignore",
                    message="Graph is not fully connected, spectral embedding" +
                    " may not work as expected.",
                    category=UserWarning)
    warnings.filterwarnings(
                    "ignore",
                    message="Array is not symmetric, and will be converted to symmetric by average with its transpose.",
                    category=UserWarning)
    warnings.filterwarnings(
                    "ignore",
                    message="Changing the sparsity structure of a csr_matrix is expensive")


    #Perform clustering
    groups = cf(**params)
    
    return groups
    



def louvain(adj_mat: np.ndarray,
                verbose: bool = False,
                shuffle_nodes: bool = False,
                tol_optimization: float = 0.0001,
                tol_aggregation: float = 0.001,
                return_adj: bool = False,
                **kwargs):
        """ Performs louvain clustering
            Returns groups, adj_mat
        """
        louvain = Louvain(modularity='dugue',
                        verbose = verbose, 
                        shuffle_nodes = shuffle_nodes, 
                        tol_optimization = tol_optimization, 
                        tol_aggregation = tol_aggregation)
        x = louvain.fit_transform(adj_mat)
        if return_adj:
            adj_mat_new = louvain.adjacency_ #Adjacency matrix between the groups
            return x, adj_mat_new
        else:
            return x
    
def spectral(adj_mat:np.ndarray,
             n: int = 200,
             n_init: int = 100,
             **kwargs):
    """Performs spectral clustering and prints to an external file. 
       If running on a pre-grouped set, include groups for correct labeling of nodes.
       Returns groups
    """

    sc = cluster.SpectralClustering(n, affinity='precomputed', n_init=n_init)
    sc.fit(adj_mat) 
    
    #Assign clusters to nodes 
    x=sc.labels_
    return x

def fluid(adj_mat: np.ndarray,
          n: int = 200,
          **kwargs):
    G = nx.from_scipy_sparse_matrix(adj_mat)
    G = G.to_undirected()
    G = nx.algorithms.community.asyn_fluidc(G, n)
    G = [x for x in G]
    x = np.empty(adj_mat.shape[0])
    for (i, ind) in enumerate(G):
        x[list(ind)] = int(i)
        
    return x

def geo_cluster(coords: np.array,
                method: str,
                n: int,
                **kwargs):
    pre = getattr(cluster,method)(n_clusters=n)
    pre.fit(coords)
    x =  pre.labels_
    return x

def geo_cluster_sc(coords: np.array,
                    adj_mat: np.ndarray,
                    method: str,
                    n: int,
                    **kwargs):
    nn = min(n, 10)
    x = geo_cluster(coords, method, nn)
    if nn < n:
        x = split_groups_sc(adj_mat,
                            groups = x, 
                            total_clusters = n, 
                            group_division = 2)
    return x

def L_L_sc(adj_mat : np.ndarray,
           n: int = 200,
           **kwargs):
    tas=np.logspace(0, -9, num=5, base=4)
    tos=np.logspace(3,-3,num=7, base=10)
    res=get_params(tas,tos,n,adj_mat)
    
    if res[0] == 'None':
        ta1=res[1]
        to1=res[2]
        ta2=None
    else:
        ta1=res[0]
        to1=1e-3
        ta2=res[1]
        to2=res[2]
    x1, adj_mat_ll = louvain(adj_mat,
                            tol_optimization=to1,
                            tol_aggregation=ta1,
                            verbose=False,
                            return_adj = True)
    adj_mat_ll.setdiag(0)
    if ta2 is not None:
        x2 =louvain(adj_mat_ll,
                                tol_optimization=to2,
                                tol_aggregation=ta2,
                                verbose=False)
        # x = np.empty(adj_mat.shape[0]).astype(int)
        # for (i, c) in enumerate(x2):
        #    x[x1 == i] = int(c)
        x = np.array([(x2[x1[i]]) for i in range(0,len(x1))]) # Convert back to node labelling
    else:
        x = x1
    x = split_groups_sc(adj_mat,
                                groups = x,
                                total_clusters=n,
                                group_division=2,
                                sort_groups=False)
    return x

def L_scn(adj_mat,
          n,
          ta = 1e-9,
          **kwargs):
    #Perform louvain with increasing ta until number of groups is
    #greater than the desired
    x1 = [0]
    while max(x1) + 1 <= n:
        x1, adj_mat_l = louvain(adj_mat,
                                    tol_optimization=1e3, #TODO should this be 1e-3?
                                    tol_aggregation=ta,
                                    verbose=False,
                                    return_adj = True)
        ta*=10
        
    # Perform spectral clustering on the resulting adjacency matrix 
    # to get the desired number of groups
    x2 = spectral(adj_mat_l, n=n)
    
    x = [(x2[x1[i]]) for i in range(0,len(x1))] # Convert back to node labelling
    # x = np.empty(adj_mat.shape[0]).astype(int)
    # for (i, c) in enumerate(x2):
    #    x[x1 == i] = int(c)
    return np.array(x)

def split_groups_sc(adj_mat:np.ndarray,
                    groups:np.ndarray,
                    total_clusters,
                    group_division:int=2,
                    sort_groups:bool=False):
    """ Performs spectral clustering within groups. Will continue to split 
        groups until the number of groups = total_clusters 
    
        total_clusters is the number of desired clusters. 
        group_division is the number of clusters into which each existing
        group will be divided
        
        Returns groups
    """
 
    while len(set(groups))<total_clusters:
        #Split the largest group
        group_num = np.argmax(np.bincount(groups))
        ind = groups == group_num
        
        # Create the adjacency matrix for group.
        adj_mat_group=adj_mat[ind][:,ind] 
        
        # perform specral clustering
        sc_groups=spectral(adj_mat_group,group_division) 
        
        #Update group numbers
        ind = np.where(ind)[0]
        for new_group in range(1,group_division):
            groups[ind[sc_groups == new_group]] = max(groups) + 1
        
    return groups


def step(tas: list, tos: list, objective: int, adj_mat:np.ndarray):
    """ Performs louvain on all ta/to combinations untill the number of
        goups is less that or equal to the objective.
        
        ta=tol_aggregation
        to=tol_optimisation
        
        Returns the combination of ta and to that produce the number of 
        groups closest to the objective (without going over)
    
    """
    ns = [] #all combos of ta,to and n
    combinations = [(x,y) for x in tas for y in tos]
    for ta, to in combinations:
        x=louvain(adj_mat,
                      tol_optimization = to,
                      tol_aggregation = ta)
        n = max(x) + 1
        ns.append(n)
    
    ns = np.array(ns)
    if min(ns) >= objective:
        ind = np.where(ns == min(ns))[0][0]
    else:
        ind = np.where(ns ==  max(ns[ns < objective]))[0][0]
        
    
    ta_min = combinations[ind][0]
    to_min = combinations[ind][1]
    n_min = ns[ind]
    return ta_min, to_min, n_min
        

def get_params(tas: list, tos: list, objective: int, adj_mat):
    """ Finds the best paramaters for louvain to get 
        number of groups<=objective. 
        Returns a tuple of (ta1, ta2 ,to2 ,n) if louvain needs to be run 
        twice, or ('None', ta1, to1, n) if only once is needed
    """
    # Perform louvain once for all combos of ta and to until the objective 
    # number of groups or less is obtained.
    res=[]
    resi=step(tas,tos,objective,adj_mat)
    resi=('None',)+resi
    #save best paramaters for single louvain
    res.append(resi)
    count=0

    while tas[count]>res[0][1] and resi[3]!=objective:
        #Perform louvain once per ta with to set to 1e-3, and then a second
        #time with all ta/to combinations
        ta=tas[count]
        n,new_adj_mat=louvain(adj_mat,
                              tol_aggregation = ta,
                              tol_optimization = 1e-3,
                              return_adj = True)
        new_adj_mat.setdiag(0)
        resi=step(tas,tos,objective,new_adj_mat)
        resi=(ta,)+resi
        #Save best paramaters for second louvain
        res.append(resi)
        count+=1
    
    # Find best paramaters overall
    ns = np.array([x[-1] for x in res])
    if min(ns) >= objective:
        ind = np.where(ns == min(ns))[0][0]
    else:
        ind = np.where(ns ==  max(ns[ns < objective]))[0][0]
    res = res[ind]

    return res