#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""This is the main code used at the MULTINETLAB for Multilayer Analysis.
   This code can create a Multilayer network object - similar to the ones in Networkx.
   It uses a Supra Adjacency Matrix as input. 
   For privacy reasons, we provide a random MST file.
"""

__author__ = "Fernando Nobrega"
__contact__ = "f.nobregasantos@amsterdamumc.nl"
__date__ = "2019/10/15"
__status__ = "Production"


####################
# Review History   #
####################

# Reviewed by Eduarda Centeno 20200909
# Reviewed by Lucas Breedt 202011..


####################
# Libraries        #
####################

####FINAL VERSIONS OF PACKAGES - PULL TO LUCAS

# Standard imports
import itertools #no attribute version
from multiprocessing import Pool #no attribute version
import time as time #no attribute version

# Third party imports
import matplotlib.pyplot as plt # version '3.1.1'
import multinetx as mx #no attribute version
import networkx as nx#version 2.3
import numpy as np#version '1.19.0'
import pandas as pd# version '0.25.1'
from pylab import pcolor, show, colorbar, xticks, yticks
import pyreadstat #version '0.2.9'
import scipy as sio  #version '1.3.1'
import scipy.io
from sklearn import preprocessing #version 0.21.3
from sklearn.preprocessing import MinMaxScaler 

### Comments with more than one '# means for further developments!!!

""" The user should define the input file (supra-adjacency matrix) in the beginning of the code. 
The code is quite robust - as long as the matrices are created using a pipeline similar to the one in the Lab."""

####################
# SETTINGS         #
####################

layer_size = 197   # Define the number of nodes per layer. We used the BNA with some regions removed

weighted = False # We are now using MST matrices. Matrices are thus not weighted - if weighted, change to True


# We now should include the file for the supra adjacency matrices here!!!

# TRAINING RANDOM MATRIX
filename = 'supra_randmst.mat'


#########################################
# CREATING LAYER TAGS                   #
#########################################


# Associating tags for each layer will be helpful for our coding. We used the ones below.
# These are the tags for the Multilayer Networks - It should match with the layers in the supra-adjacency matrix.
 
print('0 = fmri, 1 = pli delta, 2 = pli theta, 3 = pli alpha1, 4 = pli alpha2, 5 = pli beta, 6 = pli gamma, 7 = DWI .') 


##IMPROVEMENT! WE CAN INCLUDE A FUNCTION TO CHECK THE TAGS FROM LUCAS FILES 

layer_tags=['0 = fmri', '1 = pli delta', '2 = pli theta', '3 = pli alpha1', '4 = pli alpha2', '5 = pli beta', '6 = pli gamma', '7 = DWI']
just_tags=['fmri', 'pli_delta', 'pli_theta', 'pli_alpha1', 'pli_alpha2', 'pli_beta', 'pli_gamma', 'DWI'] 
plot_tags=['fMRI', 'PLI delta', 'PLI theta', 'PLI alpha1', 'PLI alpha2', 'PLI beta', 'PLI gamma', 'DWI'] 



Layer_dic = {}
for i in range(0 , len(just_tags)):
    Layer_dic[i] = just_tags[i]
print(Layer_dic)


#############################################
# LOADING THE MATRICES                      #
############################################# 

# This will be where are the real data for computing all Multilayer functions!!! 

# Notice that, from now on, every function received the data as input

# If one create other data such as random, weighed, etc, you just need to include it here !!!!

# ATTENTION!!! THIS IS THE OBJECT YOU ARE GOING TO USE FOR THE REMAINING OF THE CODE!!!!

Supra_MST = scipy.io.loadmat(filename)


##IMPROVEMENT - INCLUDE VERBOSE FUNCTION TO MAKE CHECKS IN THE CODE!!!



#######################
# SANITY CHECK        #
#######################

#Check that this is in fact a dicionary
print(type(Supra_MST))

#print(type(Supra_Weighted)) #Shows that This is in fact a dicionary

print(Supra_MST.keys())



###################################
# PREPARING THE MULTILAYER        #
###################################

###IMPROVEMENT - make it independent of the layer size!!!
  
def Prepare_Multilayer(Data,list_of_layers,N=layer_size):
    """Converts the data to a Multinetx friendly object based in the inputed data.
    Parameters
    ----------
    Data : A preloaded .mat  - Ex: Supra_MST
    
    
    list_of_layers: a list of numbers corresponding to the Multilayer you want to create - 
    Ex: If you want a Multilayer with fmri, pli_delta, pli theta, and pli beta using the tags: 0=fmri', '1=pli delta', '2= pli theta','5 = pli beta'
    list_of_layers = [0,1,2,5]
    
    Returns
    -------
    out : list of layers.    
    Note: This is convenient - since we can have a Database with all 14 layers, but we may want to study only a smaller number of layers. 
    In short, the layers are based in the supra-adjacency matrices - but you can choose exactly the layers you want here"
    """
    
    
    #In the matlab file the element [-1] gives the matrices
    name = list(Data.keys())[-1]
    multilayer = Data[name]

    # Just checking if there are NaNs
    where_are_NaNs = np.isnan(multilayer)
    multilayer[where_are_NaNs] = 0
    #layer_size = 197 # This are the numbers of nodes in the layer
    #N = layer_size
    layer_list = list_of_layers

    layers=[]
    for i in layer_list:
        layers.append(multilayer[(i*N):(i+1)*N,(i*N):(i+1)*N,:])
    return layers


# This creates a multilayer network for each individual (This is the new one)
def multlayerG(individual, Data, list_of_single_layers,N=layer_size):
    """Creates a Multilayer Network for an individual, given the data, and a list of layers.
    Parameters
    ----------
    
    Individual: an integer from [0, S-1], where S is the size of the cohort.
    
    Data : A preloaded .mat  - Ex: Supra_MST
    
    
    list_of_layers: a list of numbers corresponding to the Multilayer you want to create - 
    Ex: If you want a Multilayer with fmri, pli_delta and pli theta, using the tags: 0=fmri', '1=pli delta', '2= pli theta'
    list_of_layers = [0,1,2]
    
    Returns
    -------
    out: A Multilayer Object for a single individual in the Data.    
   
    """
    
    
    "Creates a multilayer for an individual i, knowing the number of layers, and the size of the layers"
    layers= Prepare_Multilayer(Data, list_of_single_layers)
    #N =197 # before was 205
    number_of_layers=len(list_of_single_layers)
    G=[]
    for j in range(0,len(list_of_single_layers)):
        "j is running over all layers for a fixed individual i"
        G.append(mx.from_numpy_matrix(layers[j][:,:,individual]))
    
        
    
# Define the type of interconnection between the layers

# This creates the supra adjacency matrix
    adj_block = mx.lil_matrix(np.zeros((N*number_of_layers,N*number_of_layers))) # N is the size of the layer

# Need to create generic adjacency blocks!!!!

# These are generic interconnection blocks!!!   
    for i in range(number_of_layers):
        for j in range(number_of_layers):
            if i == j:
                adj_block[i*N:  (i+1)*N,  j*N:(j+1)*N] = np.zeros(N)
            else:
                adj_block[i*N:  (i+1)*N,  j*N:(j+1)*N] = np.identity(N)    

    mg = mx.MultilayerGraph(list_of_layers=G,inter_adjacency_matrix=adj_block)
    mg.set_edges_weights(intra_layer_edges_weight=1,inter_layer_edges_weight=1)

 
    
    return mg


#############################
# CREATING THE AGGREGATE    #
#############################
    
# ATENTION: THERE ARE SEVERAL OPTIONS HERE - WE ARE USING A SIMILAR ONE DONE BY MUXVIx  at http://muxviz.net/tutorial.php  
# THIS IS AN INTERMEDIATE FUNCTION SO THAT THE OUTPUT OF THE OTHER FUNCTIONS ARE PRINTED 'PER NODE'

# This is the aggregate I suppose Muxviz is using - we double-checked this with Lucas - There are other ways of creating the aggregate.

def MVaggregate(multiple_layers_list, number_layers):
    """Creates an aggregate output from a Multilayer Network, given a multiple_layers_list and the number of layers
    Parameters
    ----------
    
    multiple_layers_list: 
    
    number_layers : 
    
    Returns
    -------
    out: An aggregate list which is the mean of the values of a Network property per node in each layer
   
    """
    k, m = divmod(len(multiple_layers_list), number_layers)
    temp = list(multiple_layers_list[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(number_layers))
    temp_mean = np.mean(temp,axis=0)
    temp_mean = temp_mean/max(temp_mean)
    #for sublists in temp:
    #    m=np.max(temp[sublists])
    #    for i in sublists:
    #        temp[sublists][i]=temp[sublists][i]/m
            
    return temp_mean




###########################################################
# HERE COMES ALL THE MULTILAYER FUNCTIONS                 #
# The strategy to create the functions are the same,      #
# We can parse all NetworkX functions here.               #
# We could try to do everything in one go, but that works #
#only if all NetworkX functions would have the same sintax#            
###########################################################

    
def Group_eigenvector_centrality(Data, list_of_single_layers):
    
    """Returns a flat list with the MVaggregate output for EC, given a Data, and a list_of_single_layers
    
    Parameters
    ----------
    Data : A preloaded .mat  - Ex: Supra_MST
    
    
    list_of_layers: a list of numbers corresponding to the Multilayer you want to create - 
    Ex: If you want a Multilayer with fmri, pli_delta, pli theta, and pli beta using the tags: 0=fmri', '1=pli delta', '2= pli theta','5 = pli beta'
    list_of_layers = [0,1,2,5]
    
    Returns
    -------
    out: An aggregate list which is the mean of the values of the EC per node in each layer
   
    """

    "This list will save all the eigenvector centralities for all individuals in all layers."
    name = list(Data.keys())[-1]

    number_of_individuals = Data[name].shape[2]

    Group_eigenvector = []
    for individual in range(number_of_individuals):
        temp = multlayerG(individual,Data,list_of_single_layers)

        m = mx.eigenvector_centrality_numpy(temp)
        #m=mx.eigenvector_centrality(multlayerG(individual,number_of_layers,list_of_layers))
        temp1 = list(m.values())
        temp2 = MVaggregate(temp1, len(list_of_single_layers)) 
        #temp2=aggregate(temp1,len(list_of_single_layers))
        # This is a list of lists with all centralities for all individuals
        Group_eigenvector.append(temp2)
        # since we want to buid a flat list 
    flat_list = [item for sublist in Group_eigenvector for item in sublist]
        
    return flat_list

def Group_clustering(Data, list_of_single_layers):
    """Returns a flat list with the aggregate output for Group clustering, given a Data, and a list_of_single_layers
    
    Parameters
    ----------
    Data : A preloaded .mat  - Ex: Supra_MST
    
    
    list_of_layers: a list of numbers corresponding to the Multilayer you want to create - 
    Ex: If you want a Multilayer with fmri, pli_delta, pli theta, and pli beta using the tags: 0=fmri', '1=pli delta', 
    '2= pli theta','5 = pli beta'
    list_of_layers = [0,1,2,5]
    
    Returns
    -------
    out: An aggregate list which is the mean of the values of the Clustering per node in each layer
   
    """
    
    name = list(Data.keys())[-1]

    number_of_individuals = Data[name].shape[2]

    Group_clustering = []
    for individual in range(number_of_individuals):
        temp = multlayerG(individual, Data, list_of_single_layers)

        m = mx.clustering(temp)
        temp1 = list(m.values())
        temp2 = MVaggregate(temp1, len(list_of_single_layers)) 
        # This is a list of lists with all centralities for all individuals
        Group_clustering.append(temp2)
    # We want to buid a flat list 
    # Check this flattened = [val for sublist in list_of_lists for val in sublist]
    flat_list = [item for sublist in Group_clustering for item in sublist]
        
    return flat_list


   
def Group_degree_centrality(Data, list_of_single_layers):
    """Returns a flat list with the aggregate output for Group degree centrality, given a Data, and a list_of_single_layers
    
    Parameters
    ----------
    Data : A preloaded .mat  - Ex: Supra_MST
    
    
    list_of_layers: a list of numbers corresponding to the Multilayer you want to create - 
    Ex: If you want a Multilayer with fmri, pli_delta, pli theta, and pli beta using the tags: 0=fmri', '1=pli delta', 
    '2= pli theta','5 = pli beta'
    list_of_layers = [0,1,2,5]
    
    Returns
    -------
    out: An aggregate list which is the mean of the values of the degree centrality per node in each layer
   
    """
    
    name = list(Data.keys())[-1]

    number_of_individuals = Data[name].shape[2]

    Group_deg_centrality = []
    for individual in range(number_of_individuals):
        temp = multlayerG(individual, Data, list_of_single_layers)

        m = mx.degree_centrality(temp)
        
        temp1 = list(m.values())
        temp2 = MVaggregate(temp1, len(list_of_single_layers)) 
        #temp2=aggregate(temp1,len(list_of_single_layers))
        # This is a list of lists with all centralities for all individuals
        Group_deg_centrality.append(temp2)
        # since we want to buid a flat list 
    # Check this flattened = [val for sublist in list_of_lists for val in sublist]
    flat_list = [item for sublist in Group_deg_centrality for item in sublist]
        
    return flat_list

    
def Group_eccentricity(Data, list_of_single_layers):
    """Returns a flat list with the aggregate output for Group eccentricity, given a Data, and a list_of_single_layers
    
    Parameters
    ----------
    Data : A preloaded .mat  - Ex: Supra_MST
    
    
    list_of_layers: a list of numbers corresponding to the Multilayer you want to create - 
    Ex: If you want a Multilayer with fmri, pli_delta, pli theta, and pli beta using the tags: 0=fmri', '1=pli delta', 
    '2= pli theta','5 = pli beta'
    list_of_layers = [0,1,2,5]
    
    Returns
    -------
    out: An aggregate list which is the mean of the values of the eccentricity per node in each layer
   
    """
    
    name = list(Data.keys())[-1]

    number_of_individuals = Data[name].shape[2]

    Group_eccentricity = []
    for individual in range(number_of_individuals):
        temp = multlayerG(individual, Data, list_of_single_layers)
        m = mx.eccentricity(temp)
        temp1 = list(m.values())
        temp2 = MVaggregate(temp1, len(list_of_single_layers)) 
        # This is a list of lists with all eccentricities for all individuals
        Group_eccentricity.append(temp2)
        # We want to buid a flat list 
    flat_list = [item for sublist in Group_eccentricity for item in sublist]
        
    return flat_list

def Non_norm_Group_eccentricity(Data, list_of_single_layers):
    """Returns a flat list with the aggregate output for Group eccentricity without normalization, given a Data, and a list_of_single_layers
    
    Parameters
    ----------
    Data : A preloaded .mat  - Ex: Supra_MST
    
    
    list_of_layers: a list of numbers corresponding to the Multilayer you want to create - 
    Ex: If you want a Multilayer with fmri, pli_delta, pli theta, and pli beta using the tags: 0=fmri', '1=pli delta', 
    '2= pli theta','5 = pli beta'
    list_of_layers = [0,1,2,5]
    
    Returns
    -------
    out: An aggregate list which is the mean of the values of the eccentricity per node in each layer without normalization
   
    """
    name = list(Data.keys())[-1]

    number_of_individuals = Data[name].shape[2]

    Group_eccentricity = []
    for individual in range(number_of_individuals):
        temp = multlayerG(individual, Data, list_of_single_layers)

        m = mx.eccentricity(temp)
        temp1 = list(m.values())
        temp2 = temp1 
        # For non normalized, we don't apply the MVaggregate function
        #MVaggregate(temp1, len(list_of_single_layers)) 
        
        #temp2=aggregate(temp1,len(list_of_single_layers))
        # This is a list of lists with all eccentricities for all individuals
        Group_eccentricity.append(temp2)
        # since we want to buid a flat list 
    flat_list = [item for sublist in Group_eccentricity for item in sublist]
        
    return flat_list

    
def Group_bet_centrality(Data, list_of_single_layers):
    """Returns a flat list with the aggregate output for Group betweeness centralities, given a Data, and a list_of_single_layers
    
    Parameters
    ----------
    Data : A preloaded .mat  - Ex: Supra_MST
    
    
    list_of_layers: a list of numbers corresponding to the Multilayer you want to create - 
    Ex: If you want a Multilayer with fmri, pli_delta, pli theta, and pli beta using the tags: 0=fmri', '1=pli delta', 
    '2= pli theta','5 = pli beta'
    list_of_layers = [0,1,2,5]
    
    Returns
    -------
    out: An aggregate list which is the mean of the values of the betweeness centralities per node in each layer
   
    """
    
    name = list(Data.keys())[-1]

    number_of_individuals = Data[name].shape[2]

    Group_bet_centrality = []
    for individual in range(number_of_individuals):
        temp = multlayerG(individual,Data,list_of_single_layers)

        m = mx.betweenness_centrality(temp)
        #m=mx.eigenvector_centrality(multlayerG(individual,number_of_layers,list_of_layers))
        temp1 = list(m.values())
        temp2 = MVaggregate(temp1, len(list_of_single_layers)) 
        #temp2=aggregate(temp1,len(list_of_single_layers))
        # This is a list of lists with all centralities for all individuals
        Group_bet_centrality.append(temp2)
        # We want to buid a flat list 
    # Check this flattened = [val for sublist in list_of_lists for val in sublist]
    flat_list = [item for sublist in Group_bet_centrality for item in sublist]
        
    return flat_list


############################################################################
#"Here comes the mean and standard deviations from all metrics we analized"#
############################################################################
    
def Group_eigenvector_centrality_mean(Data, list_of_single_layers):
    
    """Returns a flat list with the aggregate output for Group eigenvector centrality mean, given a Data, and a list_of_single_layers
    
    Parameters
    ----------
    Data : A preloaded .mat  - Ex: Supra_MST
    
    
    list_of_layers: a list of numbers corresponding to the Multilayer you want to create - 
    Ex: If you want a Multilayer with fmri, pli_delta, pli theta, and pli beta using the tags: 0=fmri', '1=pli delta', 
    '2= pli theta','5 = pli beta'
    list_of_layers = [0,1,2,5]
    
    Returns
    -------
    out: An aggregate list which is the mean of the values of the betweeness centralities per node in each layer
   
    """
    
    "This function returns the Group eigenvector centrality mean for all individuals"
    name = list(Data.keys())[-1]

    number_of_individuals = Data[name].shape[2]

    Group_eigenvector_mean = []
    for individual in range(number_of_individuals):
        temp = multlayerG(individual, Data, list_of_single_layers)

        m = mx.eigenvector_centrality_numpy(temp)
        temp1 = list(m.values())
        temp2 = MVaggregate(temp1, len(list_of_single_layers)) 
        # Now we just compute the mean
        Group_eigenvector_mean.append(np.mean(temp2))
        
    return (Group_eigenvector_mean)





def Group_eigenvector_centrality_std(Data, list_of_single_layers):
        
    """Returns a flat list with the aggregate output for Group eigenvector centrality standard deviation, given a Data, and a list_of_single_layers
    
    Parameters
    ----------
    Data : A preloaded .mat  - Ex: Supra_MST
    
    
    list_of_layers: a list of numbers corresponding to the Multilayer you want to create - 
    Ex: If you want a Multilayer with fmri, pli_delta, pli theta, and pli beta using the tags: 0=fmri', '1=pli delta', 
    '2= pli theta','5 = pli beta'
    list_of_layers = [0,1,2,5]
    
    Returns
    -------
    out: An aggregate list which is the mean of the values of the betweeness centralities per node in each layer
   
    """
    
    "This function returns the Group eigenvector centrality standard deviation for all individuals"
    
    name = list(Data.keys())[-1]

    number_of_individuals = Data[name].shape[2]
    Group_eigenvector_std = []
    for individual in range(number_of_individuals):
        temp = multlayerG(individual, Data, list_of_single_layers)
        m = mx.eigenvector_centrality_numpy(temp)
        #m=mx.eigenvector_centrality(multlayerG(individual,number_of_layers,list_of_layers))
        temp1 = list(m.values())
        temp2 = MVaggregate(temp1, len(list_of_single_layers)) 
        # This is MV aggregate - we can change then later for something else if needed
        
        Group_eigenvector_std.append(np.std(temp2))
        
    return (Group_eigenvector_std)

def Eigenvector_centrality(individual, Data, list_of_single_layers):
    """Returns a histogram with the values of the Eigenvector centrality for all nodes for a chosen individual."
    Parameters
    ----------
    individual: an integer from [0, S-1], where S is the size of the cohort.
        
    Data : A preloaded .mat  - Ex: Supra_MST
    
    
    list_of_layers: a list of numbers corresponding to the Multilayer you want to create - 
    Ex: If you want a Multilayer with fmri, pli_delta and pli theta, using the tags: 0=fmri', '1=pli delta', '2= pli theta'
    list_of_layers = [0,1,2]
    
    Returns
    -------
    out: A list with EC for one individual
   
    """
    print('layers =',[layer_tags[i] for i in list_of_single_layers])
    m = mx.eigenvector_centrality_numpy(multlayerG(individual, Data, list_of_single_layers))
    temp1 = list(m.values())
    temp2=MVaggregate(temp1,len(list_of_single_layers))
   
    return temp2

def Group_degree_centrality_mean(Data,list_of_single_layers):
     """Returns a flat list with the aggregate output for Group degree centrality mean, given a Data, and a list_of_single_layers
    
    Parameters
    ----------
    Data : A preloaded .mat  - Ex: Supra_MST
    
    
    list_of_layers: a list of numbers corresponding to the Multilayer you want to create - 
    Ex: If you want a Multilayer with fmri, pli_delta, pli theta, and pli beta using the tags: 0=fmri', '1=pli delta', 
    '2= pli theta','5 = pli beta'
    list_of_single_layers  = [0,1,2,5]
    
    Returns
    -------
    out: An aggregate list which is the mean of the values of the betweeness centralities per node in each layer
   
    """
    
#    "This function returns the Group degree centrality mean for all individuals"
    
     print('layers =',[layer_tags[i] for i in list_of_single_layers])
     name = list(Data.keys())[-1]
     number_of_individuals = Data[name].shape[2]

     "This list will save all the eigenvector centralities means for all individuals"
     Group_degree_centrality_mean = []
     for individual in range(number_of_individuals):
         m = mx.degree_centrality(multlayerG(individual,Data,list_of_single_layers))
         temp1 = list(m.values()) # this is not aggregated 
         temp2 = MVaggregate(temp1, len(list_of_single_layers)) # This is Mux Viz aggregate
        #temp2=aggregate(temp1,len(list_of_single_layers))
        # IF you want - by any chance to do a different agreggate you should change the line above
         Group_degree_centrality_mean.append(np.mean(temp2))
         
     return (Group_degree_centrality_mean)



def Group_degree_centrality_std(Data, list_of_single_layers):
    """Returns a flat list with the aggregate output for Group degree centrality std, given a Data, and a list_of_single_layers
    
    Parameters
    ----------
    Data : A preloaded .mat  - Ex: Supra_MST
    
    
    list_of_layers: a list of numbers corresponding to the Multilayer you want to create - 
    Ex: If you want a Multilayer with fmri, pli_delta, pli theta, and pli beta using the tags: 0=fmri', '1=pli delta', 
    '2= pli theta','5 = pli beta'
    list_of_single_layers = [0,1,2,5]
    
    Returns
    -------
    out: An aggregate list which is the mean of the values of the betweeness centralities per node in each layer
   
    """
    print('layers =',[layer_tags[i] for i in list_of_single_layers])
    name = list(Data.keys())[-1]
    number_of_individuals = Data[name].shape[2]
    "This list will save all the eigenvector centralities means for all individuals"
    Group_degree_centrality_std = []
    for individual in range(number_of_individuals):
        m = mx.degree_centrality(multlayerG(individual, Data, list_of_single_layers))
        temp1 = list(m.values())
        temp2 = MVaggregate(temp1, len(list_of_single_layers))  # This is Mux Viz aggregate
        #temp2=aggregate(temp1,len(list_of_single_layers)) You can change the aggregate here
        print(temp2)
        Group_degree_centrality_std.append(np.std(temp2))
        
    return (Group_degree_centrality_std)


###################################################
# SOME PLOTTING FUNCTIONS FOR MULTILAYER ANALYSIS #  
###################################################



def Plot_Group_EC(Data, list_of_single_layers):
    
    """Returns a histogram plot with the values of the Eigenvalue centrality for all nodes across all individuals."
    Parameters
    ----------
        
    Data : A preloaded .mat  - Ex: Supra_MST
    
    
    list_of_layers: a list of numbers corresponding to the Multilayer you want to create - 
    Ex: If you want a Multilayer with fmri, pli_delta and pli theta, using the tags: 0=fmri', '1=pli delta', '2= pli theta'
    list_of_layers = [0,1,2]
    
    Returns
    -------
    out : A histogram plot with EC for all nodes across all individuals
   
    """
    
    
    
    print('layers =',[layer_tags[i] for i in list_of_single_layers])

    temp = Group_eigenvector_centrality(Data,list_of_single_layers)
    plt.figure(figsize=(8,5))
    plt.hist(temp)
    # We can edit here the output if we have a vector with the name of the layers
    plt.xlabel('Eig. centr. - aggr- all nodes all individuals ', fontsize=20)
    #plt.xlim(-5,220)
    plt.ylabel("frequence", fontsize=20)
    #plt.xlim(40, 160)
    plt.ylim(0, 3500)
    #plt.title('individual '+str(individual))
    plt.show()


def Plot_EC(individual, Data, list_of_single_layers):
    """Returns a histogram with the values of the Eigenvalue centrality for all nodes for a chosen individual."
    Parameters
    ----------
    individual: an integer from [0, S-1], where S is the size of the cohort.
        
    Data : A preloaded .mat  - Ex: Supra_MST
    
    
    list_of_layers: a list of numbers corresponding to the Multilayer you want to create - 
    Ex: If you want a Multilayer with fmri, pli_delta and pli theta, using the tags: 0=fmri', '1=pli delta', '2= pli theta'
    list_of_layers = [0,1,2]
    
    Returns
    -------
    out: A histogram with EC for one individual
   
    """
    
    print('layers =',[layer_tags[i] for i in list_of_single_layers])
    #multlayerG(individual,Data,list_of_single_layers)
    m = mx.eigenvector_centrality_numpy(multlayerG(individual,Data,list_of_single_layers))
    #temp=multlayer3(i)
    temp1 = list(m.values())
    temp2 = MVaggregate(temp1, len(list_of_single_layers)) 
    # This is the Mux Viz aggregate - We change the aggregate here if yo want later
    plt.hist((temp2))
    ###IMPROVEMENT: We can edit here the output if we have a vector with the name of the layers
    plt.xlabel('Eigenvector centrality - aggregate ', fontsize=20)
    #We may also want to choose the range for the plot.
    #plt.xlim(min,max)
    plt.ylabel("frequence", fontsize=20)
    plt.ylim(0, 100)
    plt.title('individual '+ str(individual))
    plt.show()
    return
             

#### IMPROVEMENT ALSO WITHIM IN THE FUTURE
def Mask_subnetwork(result, target):
    """Returns a multilayer metric narrowed for a given list of nodes, which for our purposes are subnetworks""
    
    Parameters
    ----------
    result: A list with the results (output) of any of Multilayer functions in this code
    
    target : A list of target nodes of interest, e.g., nodes from DFN or FPN
    
    Returns
    -------
    out: A list for the results narrowed for the target nodes, i.e., If you say the target nodes for a given subnetwork, this function returns only the results of the list associated with the target nodes"
    """

    chunks = [result[x:x+layer_size] for x in range(0, len(result), layer_size)]
    mask = [chunk[x] for chunk in chunks for x in target]
    
    return mask

# This creates a SPSS file - Data is the data, name is the name of the file and tag is the collum name
def SaveSPSS(Data, name, tag):
    """Returns a .csv file for further analysis using SPSS
    Parameters
    ----------
    Data: The desired Data you want to save
    name: The name of the file you want to save
    tag: The tag for the variable/column in your .csv file
    """
    # Obs: Notice that if you want to get results only for a subnetwork, we should first do:
    #Data=Mask_subnetwork(result,target)
    #before saving this file
    cols = [tag]
    df = pd.DataFrame(Data, columns=cols)
    df.to_csv(name+'.csv')
    
    return


# Check if You wanna put Mask true of False - create a bollean that does stuff when a mask is choosen or not
def Function_output(function, Data, filename, colname, layers,N=layer_size):
    """Returns the desired output for the MumoBrain database, or any other database organized similarly
    
    THIS IS PROBLABLY THE MOST IMPORTANT FUNCTION FOR THE USER OF THIS CODE, SINCE EVERYTHING WAS BUILT TO REACH THIS STAGE HERE
    
    Parameters
    ----------
    function: One of the functions developed in this code for Multilayer Networks
    Data: The Data we want to use: e.g., Supra_MST
    filename: The name of the file you want to save
    colname: the name of the column tag in your file
    layers: list of desired layers.

    """
    temp = function(Data, layers)
    #You should include the desired subnetwork here. Now we have the whole Network.
    Sub_Net = list(range(0,N)) 
    # Ex: If you want FPN, Sub_Net=[16,17,18,19,20,21,28,29,30,31,93,94,123,124,133,134,163,164]
    
    temp_Sub_Net = Mask_subnetwork(temp, Sub_Net)
    SaveSPSS(temp_Sub_Net, filename, colname)
    
    return




