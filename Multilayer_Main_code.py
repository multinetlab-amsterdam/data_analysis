
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 15 11:31:56 2019

@author: Turing
"""

"importing important stuff."

""
import numpy as np
import multinetx as mx
import networkx as nx
import matplotlib.pyplot as plt

import scipy as sio 
import scipy.io
from sklearn import preprocessing
from pylab import pcolor, show, colorbar, xticks, yticks
import networkx as nx
import matplotlib.pyplot as plt
import itertools
from sklearn.preprocessing import MinMaxScaler
import time as time
from multiprocessing import Pool
import pyreadstat
import pandas as pd

"This main code used at MULTINETLAB for Multilayer Analysis - This code can create a Multilayer network object - similar to the ones in Networkx - having as an input a Supra Adjacency Matrix."

"For privacy reasons, we used random MST matrices on Github " 


# One should just change the value of (which is now  N=203) as the size of the matrix according to the demand

"The user should declare the supra matrices files in the beginning of the code - The code is quite robust - as long as the matrices are created using a similar pipeline as the one in the Lab."

####################
# SOME SETTINGS!!!!!
####################

layer_size=197   # THis was our setting for BNA ATLAS - for AAL we did 75!

weighted=False # We are using MST's so the matrices are not weighted - if the matrices are weighted, we should say True!!!

#Should check which file I should put here!!!!
# We now should include the file for the supra adjacency matrixed here!!!

# TRAINING RANDOM MATRIX

filename1='supra_randmst.mat'

# SOME NAMES USED IN THE LAB BEFORE!!!
#filename1='aal_supra_mst_full.mat'#'supra_MST_v2'#'supra_MST.mat' #14 layers including dwi
#filename2='aal_rand_supra_mst.mat'#'supra_MST_randomized.mat'
#filename3='supra_weighted.mat'
#filename1='supra_MST_full_203rois.mat'
#filename1='supra_MST_first25.mat'
#filename2='supra_first25.mat'
#filename4='supra_weighted_randomized.mat'


#########################################
#CREATING TAGS FOR THE LAYERS
#########################################



# Associating tags for each layer will be helpful for our coding. We used one's bellow
# These are the tags for the Multilayer Networks - It should match the with the layers in the SupraAdjacency Matrix. 
print('0 = fmri, 1=pli delta, 2= pli theta, 3= pli alpha1, 4= pli alpha2, 5 = pli beta, 6 = pli gamma, 7 = DWI.') 


#IMPROVEMENT!!! INCLUDE A FUNCTION TO CHECK THE TAGS FROM LUCAS FILES 

layer_tags=['0=fmri', '1=pli delta', '2= pli theta', '3= pli alpha1', '4= pli alpha2', '5 = pli beta', '6 = pli gamma', '7 = DWI']
just_tags=['fmri', 'pli_delta', 'pli_theta', 'pli_alpha1', 'pli_alpha2', 'pli_beta', 'pli_gamma', 'DWI'] 
plot_tags=['fMRI', 'PLI delta', 'PLI theta', 'PLI alpha1', 'PLI alpha2', 'PLI beta', 'PLI gamma', 'DWI'] 



Layer_dic={}
for i in range(0,len(just_tags)):
    Layer_dic[i]=just_tags[i]
print(Layer_dic)


#############################################
#LOADING THE MATRICES
############################################# 

   

#This is the real Data for all Multilayer functions!!! 
# Notice that, from now on, every function received the data as input
#If one create other data such as random, weighed, etc, you just need to include it here !!!!

#THIS IS THE OBJECT YOU ARE GOING TO USE FOR THE REMMAINING OF THE CODE!!!!

Supra_MST = scipy.io.loadmat(filename1)

## OTHER DATA GENERATED ARE WRITTEN BELLOW

#Supra_MST_random =scipy.io.loadmat(filename2)
#Supra_weighted =scipy.io.loadmat(filename3)
#Supra_weighted_random =scipy.io.loadmat(filename4)

###IMPROVEMENT - INCLUDE VERBOSE TO MAKE CHECKS IN THE CODE!!!



#######################
#SANITY CHECK
#######################

#Check that this is in fact a dicionary
print(type(Supra_MST))

#print(type(Supra_Weighted)) #Shows that This is in fact a dicionary

print(Supra_MST.keys())

# The command bellow gives the keys
#print(Supra_Weighted.keys())# The command bellow gives the keys

#multilayer=mat['aal_multilayer_adjmats']

#def Create_Multilayer(Data):
  
###IMPROVEMENT - make it independent of the layer size!!!
  
def Prepare_Multilayer(Data,list_of_layers):
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
    
    
    #In the matlab file the element [-1] gives the matrices#
    name=list(Data.keys())[-1]
    multilayer=Data[name]
    # Just checking if there are nan
        
    where_are_NaNs = np.isnan(multilayer)
    multilayer[where_are_NaNs] = 0
    layer_size=197#197 # This are the numbers of nodes in the layer
    N=layer_size
    layer_list=list_of_layers

    layers=[]
    for i in layer_list:
        layers.append(multilayer[(i*N):(i+1)*N,(i*N):(i+1)*N,:])
    return layers





"This creates a multilayer network for each individual"
# This is the new one
def multlayerG(individual,Data,list_of_single_layers):
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
    layers= Prepare_Multilayer(Data,list_of_single_layers)
    N =197# 197 # before was 205
    number_of_layers=len(list_of_single_layers)
    G=[]
    for j in range(0,len(list_of_single_layers)):
        "j is running over all layers for a fixed individual i"
        G.append(mx.from_numpy_matrix(layers[j][:,:,individual]))
    
        
    

#Define the type of interconnection between the layers
#"N is the size of the layer""

# This creates the supra adjacency matrix"
    adj_block = mx.lil_matrix(np.zeros((N*number_of_layers,N*number_of_layers)))

#Need to create generic adjacency blocks!!!!

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
#CREATING THE AGGREGATE
#############################
    
#ATENTION: THERE ARE SEVERAL OPTIONS HERE - WE ARE USING A SIMILAR ONE DONE BY MUXVIS - THIS WILL BE AN INTERMEDIATE FUNCTION SO THAT THE OUTPUT OF THE OTHER FUNCTIONS COMES BACK PER NODE

# This is the aggregate I suppose Muxviz is using - we double-checked this (with Lucas)

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
    temp=list(multiple_layers_list[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(number_layers))
    temp_mean=np.mean(temp,axis=0)
    temp_mean=temp_mean/max(temp_mean)
    #for sublists in temp:
    #    m=np.max(temp[sublists])
    #    for i in sublists:
    #        temp[sublists][i]=temp[sublists][i]/m
            
    return temp_mean

#def newaggregate(multiple_layers_list, number_layers):
#    k, m = divmod(len(multiple_layers_list), number_layers)
#    temp=list(multiple_layers_list[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(number_layers))
    #for sublists in temp:
        #m=np.max(temp[sublists])
        #for i in sublists:
            #temp[sublists][i]=temp[sublists][i]/m
            
#    return temp#np.mean(temp, axis=0)

    
#def aggregate(multiple_layers_list, number_layers):
#    k, m = divmod(len(multiple_layers_list), number_layers)
#    temp=list(multiple_layers_list[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(number_layers))
#    return np.mean(temp, axis=0)



#def multlayerG(individual,number_of_layers,list_of_layers):
    

    
# Now we create the group eigenvector centrality for the whole group  
    

# This function is useful for the aggregate - It is the missing step between networkx and Multilayer


"Here come all the functions - the strategy to create the functions are the same - we can parse all NetworkX functions here."



#########################################
#HERE COMES ALL THE MULTILAYER FUNCTIONS
#########################################



    
def Group_eigenvector_centrality(Data,list_of_single_layers):
    
    """Returns a flat list with the aggregate output for EC, given a Data, and a list_of_single_layers
    
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
    name=list(Data.keys())[-1]

    number_of_individuals=Data[name].shape[2]

    Group_eigenvector=[]
    for individual in range(number_of_individuals):
        temp=multlayerG(individual,Data,list_of_single_layers)

        m=mx.eigenvector_centrality_numpy(temp)
        #m=mx.eigenvector_centrality(multlayerG(individual,number_of_layers,list_of_layers))
    #temp=multlayer3(i)
        temp1=list(m.values())
        temp2=MVaggregate(temp1, len(list_of_single_layers)) 
        #temp2=aggregate(temp1,len(list_of_single_layers))
        # This is a list of lists with all centralities for all individuals
        Group_eigenvector.append(temp2)
        # since we want to buid a flat list 
    flat_list = [item for sublist in Group_eigenvector for item in sublist]
        
    return flat_list

def Group_clustering(Data,list_of_single_layers):
    
    """Returns a flat list with the aggregate output for Group clustering, given a Data, and a list_of_single_layers
    
    Parameters
    ----------
    Data : A preloaded .mat  - Ex: Supra_MST
    
    
    list_of_layers: a list of numbers corresponding to the Multilayer you want to create - 
    Ex: If you want a Multilayer with fmri, pli_delta, pli theta, and pli beta using the tags: 0=fmri', '1=pli delta', '2= pli theta','5 = pli beta'
    list_of_layers = [0,1,2,5]
    
    Returns
    -------
    out: An aggregate list which is the mean of the values of the Clustering per node in each layer"""

    
    name=list(Data.keys())[-1]
    number_of_individuals=Data[name].shape[2]

    Group_clustering=[]
    for individual in range(number_of_individuals):
        temp=multlayerG(individual,Data,list_of_single_layers)

        m=mx.clustering(temp)
        #m=mx.eigenvector_centrality(multlayerG(individual,number_of_layers,list_of_layers))
    #temp=multlayer3(i)
        temp1=list(m.values())
        temp2=MVaggregate(temp1, len(list_of_single_layers)) 
        #temp2=aggregate(temp1,len(list_of_single_layers))
        # This is a list of lists with all centralities for all individuals
        Group_clustering.append(temp2)
        # since we want to buid a flat list 
    # Check this flattened = [val for sublist in list_of_lists for val in sublist]
    flat_list = [item for sublist in Group_clustering for item in sublist]
        
    return flat_list


   
def Group_degree_centrality(Data,list_of_single_layers):
    "This list will save all the betweeness centralities means for all individuals"
    name=list(Data.keys())[-1]

    number_of_individuals=Data[name].shape[2]

    Group_deg_centrality=[]
    for individual in range(number_of_individuals):
        temp=multlayerG(individual,Data,list_of_single_layers)

        m=mx.degree_centrality(temp)
        #m=mx.eigenvector_centrality(multlayerG(individual,number_of_layers,list_of_layers))
    #temp=multlayer3(i)
        temp1=list(m.values())
        temp2=MVaggregate(temp1, len(list_of_single_layers)) 
        #temp2=aggregate(temp1,len(list_of_single_layers))
        # This is a list of lists with all centralities for all individuals
        Group_deg_centrality.append(temp2)
        # since we want to buid a flat list 
    # Check this flattened = [val for sublist in list_of_lists for val in sublist]
    flat_list = [item for sublist in Group_deg_centrality for item in sublist]
        
    return flat_list

    
def Group_eccentricity(Data,list_of_single_layers):
    #m=multlayerGNew(3,Supra_MST,[0,5,1])
    "This list will save all the eccentricity for all individuals"
    name=list(Data.keys())[-1]

    number_of_individuals=Data[name].shape[2]

    Group_eccentricity=[]
    for individual in range(number_of_individuals):
        temp=multlayerG(individual,Data,list_of_single_layers)

        m=mx.eccentricity(temp)
        #m=mx.eigenvector_centrality(multlayerG(individual,number_of_layers,list_of_layers))
    #temp=multlayer3(i)
        temp1=list(m.values())
        temp2=MVaggregate(temp1, len(list_of_single_layers)) 
        #temp2=aggregate(temp1,len(list_of_single_layers))
        # This is a list of lists with all centralities for all individuals
        Group_eccentricity.append(temp2)
        # since we want to buid a flat list 
    flat_list = [item for sublist in Group_eccentricity for item in sublist]
        
    return flat_list

def Non_norm_Group_eccentricity(Data,list_of_single_layers):
    #m=multlayerGNew(3,Supra_MST,[0,5,1])
    "This list will save all the nonormalized eccentricity for all individuals"
    name=list(Data.keys())[-1]

    number_of_individuals=Data[name].shape[2]

    Group_eccentricity=[]
    for individual in range(number_of_individuals):
        temp=multlayerG(individual,Data,list_of_single_layers)

        m=mx.eccentricity(temp)
        #m=mx.eigenvector_centrality(multlayerG(individual,number_of_layers,list_of_layers))
    #temp=multlayer3(i)
        temp1=list(m.values())
        temp2=temp1#MVaggregate(temp1, len(list_of_single_layers)) 
        
        #temp2=aggregate(temp1,len(list_of_single_layers))
        # This is a list of lists with all centralities for all individuals
        Group_eccentricity.append(temp2)
        # since we want to buid a flat list 
    flat_list = [item for sublist in Group_eccentricity for item in sublist]
        
    return flat_list

    
def Group_bet_centrality(Data,list_of_single_layers):
    "This list will save all the betweeness centralities means for all individuals"
    name=list(Data.keys())[-1]

    number_of_individuals=Data[name].shape[2]

    Group_bet_centrality=[]
    for individual in range(number_of_individuals):
        temp=multlayerG(individual,Data,list_of_single_layers)

        m=mx.betweenness_centrality(temp)
        #m=mx.eigenvector_centrality(multlayerG(individual,number_of_layers,list_of_layers))
    #temp=multlayer3(i)
        temp1=list(m.values())
        temp2=MVaggregate(temp1, len(list_of_single_layers)) 
        #temp2=aggregate(temp1,len(list_of_single_layers))
        # This is a list of lists with all centralities for all individuals
        Group_bet_centrality.append(temp2)
        # since we want to buid a flat list 
    # Check this flattened = [val for sublist in list_of_lists for val in sublist]
    flat_list = [item for sublist in Group_bet_centrality for item in sublist]
        
    return flat_list



"Here comes the mean and standard deviations from all metrics we analized"
    
def Group_eigenvector_centrality_mean(Data,list_of_single_layers):
    #m=multlayerGNew(3,Supra_MST,[0,5,1])
    "This list will save all the eigenvector centralities means for all individuals"
    name=list(Data.keys())[-1]

    number_of_individuals=Data[name].shape[2]

    Group_eigenvector_mean=[]
    for individual in range(number_of_individuals):
        temp=multlayerG(individual,Data,list_of_single_layers)

        m=mx.eigenvector_centrality_numpy(temp)
        #m=mx.eigenvector_centrality(multlayerG(individual,number_of_layers,list_of_layers))
    #temp=multlayer3(i)
        temp1=list(m.values())
        temp2=MVaggregate(temp1, len(list_of_single_layers)) 
        #temp2=aggregate(temp1,len(list_of_single_layers))
        
        Group_eigenvector_mean.append(np.mean(temp2))
        
    return (Group_eigenvector_mean)





def Group_eigenvector_centrality_std(Data,list_of_single_layers):
    "This list will save all the eigenvector centralities stds for all individuals"
    name=list(Data.keys())[-1]

    number_of_individuals=Data[name].shape[2]
    Group_eigenvector_std=[]
    for individual in range(number_of_individuals):
        temp=multlayerG(individual,Data,list_of_single_layers)
        m=mx.eigenvector_centrality_numpy(temp)
        #m=mx.eigenvector_centrality(multlayerG(individual,number_of_layers,list_of_layers))
    #temp=multlayer3(i)
        temp1=list(m.values())
        temp2=MVaggregate(temp1, len(list_of_single_layers)) # This is MV aggregate - we can change then later for something else
        #temp2=aggregate(temp1,len(list_of_single_layers))
        
        Group_eigenvector_std.append(np.std(temp2))
        
    return (Group_eigenvector_std)


###############################
#SOME PLOTTING FUNCTIONS FOR EC   
###############################
#"Those are functions to plot the previous functions."


def Plot_Group_EC(Data,list_of_single_layers):
    """This function plots a histogram with the values of the Eigenvalue centrality for all nodes across all individuals."

    Parameters
    ----------
        
    Data : A preloaded .mat  - Ex: Supra_MST
    
    
    list_of_layers: a list of numbers corresponding to the Multilayer you want to create - 
    Ex: If you want a Multilayer with fmri, pli_delta and pli theta, using the tags: 0=fmri', '1=pli delta', '2= pli theta'
    list_of_layers = [0,1,2]
    
    Returns
    -------
    out : A histogram with EC for all nodes across all individuals
   
    """
    
    
    
    print('layers =',[layer_tags[i] for i in list_of_single_layers])

    temp=Group_eigenvector_centrality(Data,list_of_single_layers)
    plt.figure(figsize=(8,5))
    plt.hist(temp)
    # We can edit here the output if we have a vector with the name of the layers
    plt.xlabel('Eig. centr. - aggr- all nodes all individuals ',fontsize=20)
    #plt.xlim(-5,220)
    plt.ylabel("frequence",fontsize=20)
    #plt.xlim(40, 160)
    plt.ylim(0, 3500)
#plt.title('individual '+str(individual))
    plt.show()


def Plot_EC(individual,Data,list_of_single_layers):
    """This function plots a histogram with the values of the Eigenvalue centrality for all nodes one individual."

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
    m=mx.eigenvector_centrality_numpy(multlayerG(individual,Data,list_of_single_layers))
    #temp=multlayer3(i)
    temp1=list(m.values())
    temp2=MVaggregate(temp1, len(list_of_single_layers)) 
    # This is the Mux Viz aggregate - We change the aggregate here if yo want later
    #temp2=aggregate(temp1,len(list_of_single_layers))
    plt.hist((temp2))
    ###IMPROVEMENT: We can edit here the output if we have a vector with the name of the layers
    plt.xlabel('Eigenvector centrality - aggregate ',fontsize=20)
    #plt.xlim(-5,220)
    plt.ylabel("frequence",fontsize=20)
    plt.ylim(0, 100)
    plt.title('individual '+str(individual))
    plt.show()
    return
    



            

#
#def aggregate(multiple_layers_list, number_layers):
#    k, m = divmod(len(multiple_layers_list), number_layers)
#    temp=list(multiple_layers_list[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(number_layers))
#    return np.mean(temp, axis=0)

    
#def aggregate(multiple_layers_list, number_layers):
    #k, m = divmod(len(multiple_layers_list), number_layers)
    #temp=list(multiple_layers_list[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(number_layers))
    #return np.mean(temp, axis=0)



def eigenvectorcentrality(individual,Data,list_of_single_layers):
    print('layers =',[layer_tags[i] for i in list_of_single_layers])
    #multlayerG(individual,Data,list_of_single_layers)
    m=mx.eigenvector_centrality_numpy(multlayerG(individual,Data,list_of_single_layers))
    #temp=multlayer3(i)
    temp1=list(m.values())
    #temp2=aggregate(temp1,len(list_of_single_layers))
    #plt.hist((temp2))
    # We can edit here the output if we have a vector with the name of the layers
    #plt.xlabel('Eigenvector centrality - aggregate ',fontsize=20)
    #plt.xlim(-5,220)
    #plt.ylabel("frequence",fontsize=20)
    #plt.title('individual '+str(individual))
    #plt.show()
    return temp1


def Group_degree_centrality_mean(Data,list_of_single_layers):
    print('layers =',[layer_tags[i] for i in list_of_single_layers])

    #print('layers =',[layer_tags[i] for i in list_of_single_layers])
    name=list(Data.keys())[-1]

    number_of_individuals=Data[name].shape[2]
    "This list will save all the eigenvector centralities means for all individuals"
    Group_degree_centrality_mean=[]
    for individual in range(number_of_individuals):
        m=mx.degree_centrality(multlayerG(individual,Data,list_of_single_layers))
    #temp=multlayer3(i)
        temp1=list(m.values()) # this is not aggregated 
        temp2=MVaggregate(temp1, len(list_of_single_layers)) # This is Mux Viz aggregate
        #temp2=aggregate(temp1,len(list_of_single_layers))
        # IF you want - by any chance to do a different agreggate you should change the line above
        Group_degree_centrality_mean.append(np.mean(temp2))
        
    return (Group_degree_centrality_mean)

def Group_degree_centrality_std(Data,list_of_single_layers):
    print('layers =',[layer_tags[i] for i in list_of_single_layers])
    name=list(Data.keys())[-1]


    number_of_individuals=Data[name].shape[2]
    "This list will save all the eigenvector centralities means for all individuals"
    Group_degree_centrality_std=[]
    for individual in range(number_of_individuals):
        m=mx.degree_centrality(multlayerG(individual,Data,list_of_single_layers))
    #temp=multlayer3(i)
        temp1=list(m.values())
        temp2=MVaggregate(temp1, len(list_of_single_layers))  # This is Mux Viz aggregate
        #temp2=aggregate(temp1,len(list_of_single_layers)) You can change the aggregate here
        print(temp2)
        Group_degree_centrality_std.append(np.std(temp2))
        
    return (Group_degree_centrality_std)
" We need to create 90 metrics - 45 for the Multilayer and 45 for the random networks."


#### IMPROVEMENT ALSO WITHIM
def Mask_subnetwork(result,target):
    "If you say the target nodes for a given subnetwork, this command you return only the results of list associated with the target nodes"
    chunks = [result[x:x+layer_size] for x in range(0, len(result), layer_size)]
    mask=[chunk[x] for chunk in chunks for x in target]
    
    return mask
# This creates a SPSS file - Data is the data, name is the name of the file and tag is the collum name
def SaveSPSS(Data,name,tag):
    #Data=Mask_subnetwork(result,target)
    cols = [tag]
    df = pd.DataFrame(Data,columns=cols)
    df.to_csv(name+'.csv')
    #pyreadstat.write_sav(df, name+'.csv')
    return


# Check if You wanna put Mask true of False - create a bollean that does stuff when a mask is choosen or not
def Function_output(function,Data,filename,colname,layers):
    # THis can be nany function we develeloped!
    temp=function(Data,layers)
    FPN=list(range(0,197)) #THIS IS WITHOUT THE MASK - before 197
    #FPN=[16,17,18,19,20,21,28,29,30,31,93,94,123,124,133,134,163,164]
    #old one - FPN=[16,17,18,19,20,21,28,29,30,31,93,94,129,130,139,140,169,170]
    temp_FPN=Mask_subnetwork(temp,FPN)
    SaveSPSS(temp_FPN,filename,colname)
    
    return




##########################################################
# HERE ARE NEARLY ALL THE RUNNING ANALYSIS I DID FOR LUCAS
##########################################################

# Create a code - called Running Stuff - So that I can just run things based in this
#Monolayer EC---------------- # Without Mask! Check function above!
#for i in range(0,8):
#    filename='Random_No_Mask_Group_MST_Mono_layer_real_'+Layer_dic[i]+'_tag_'+str(i)
#    colname='Rando,_No_Mask_Group_MST_Mono_layer_real_'+Layer_dic[i]+'_tag_'+str(i)
#    print(filename)
#    print(colname)
#    function=Group_eigenvector_centrality
#    Data=Supra_MST
#    Function_output(function,Data,filename,colname,[i])
#    print('we did it')






#M0no RANDOM_EC-----------------
#for i in range(0,8):
#    filename='EC_No_Mask_Group_MST_Mono_layer_random_'+Layer_dic[i]+'_tag_'+str(i)
#    colname='EC_No_Mask_Group_MST_Mono_layer_random_'+Layer_dic[i]+'_tag_'+str(i)
#    print(filename)
#    print(colname)
#    function=Group_eigenvector_centrality
#    Data=Supra_MST_random
#    Function_output(function,Data,filename,colname,[i])
#    print('we did it')
    

#
#MULTILAYER _EC-----real ------
#filename='EC_No_Mask_Group_MST_Multi_layer_real_'#+Layer_dic[i]+'_tag_'+str(i)
#colname='EC_No_Mask_Group_MST_Multi_layer_real_'#+Layer_dic[i]+'_tag_'+str(i)
#print(filename)
#print(colname)
#function=Group_eigenvector_centrality
#Data=Supra_MST
#Function_output(function,Data,filename,colname,list(range(8)))
#print('we did it')

#MULTI_EC------random

#filename='EC_No_Mask_Group_MST_Multi_layer_random_'#+Layer_dic[i]+'_tag_'+str(i)
#colname='EC_No_Mask_Group_MST_Multi_layer_random_'#+Layer_dic[i]+'_tag_'+str(i)
#print(filename)
#print(colname)
#function=Group_eigenvector_centrality
#Data=Supra_MST_random
#Function_output(function,Data,filename,colname,list(range(8)))
#print('we did it')
    
#___________________________________________________________
    #-------------------
    
#MULTILAYER clustering
#filename='Clustering_no_mask_Group_MST_Multi_layer_real_'#+Layer_dic[i]+'_tag_'+str(i)
#colname='Clustering_no_mask_Group_MST_Multi_layer_real_'#+Layer_dic[i]+'_tag_'+str(i)
#print(filename)
#print(colname)
#function=Group_clustering
#Data=Supra_MST
#Function_output(function,Data,filename,colname,list(range(8)))
#print('we did it')
#------------------------------------------------------------
    
    

#MULTILAYER Degree Centrality
#filename='Degree_cent_no_mask_Group_MST_Multi_layer_real_'#+Layer_dic[i]+'_tag_'+str(i)
#colname='Degree_cent_no_mask_Group_MST_Multi_layer_real_'#+Layer_dic[i]+'_tag_'+str(i)
#print(filename)
#print(colname)
#function=Group_degree_centrality
#Data=Supra_MST
#Function_output(function,Data,filename,colname,list(range(8)))
#print('we did it')
#------------------------------------------------------------


   
#mono Group_bet_centrality

#for i in range(0,8):
#    filename='Bet_cent_FPN_Group_MST_Mono_layer_real_'+Layer_dic[i]+'_tag_'+str(i)
#    colname='Bet_cent_FPN_Group_MST_Mono_layer_real_'+Layer_dic[i]+'_tag_'+str(i)
#    print(filename)
#    print(colname)
#    function=Group_bet_centrality
#    Data=Supra_MST
#    Function_output(function,Data,filename,colname,[i])
#    print('we did it')
#----------------------------

# mono BET_CENT_RANDOM
#for i in range(0,8):
#    filename='Bet_cent_FPN_Group_MST_Mono_layer_random_'+Layer_dic[i]+'_tag_'+str(i)
#    colname='Bet_cent_FPN_Group_MST_Mono_layer_random_'+Layer_dic[i]+'_tag_'+str(i)
#    print(filename)
#    print(colname)
#    function=Group_bet_centrality
#    Data=Supra_MST_random
#    Function_output(function,Data,filename,colname,[i])
#    print('we did it')



#-------------------
#MULTILAYER _BC
#filename='Bet_cent_FPN_Group_MST_Multi_layer_real_'#+Layer_dic[i]+'_tag_'+str(i)
#colname='Bet_cent_FPN_Group_MST_Multi_layer_real_'#+Layer_dic[i]+'_tag_'+str(i)
#print(filename)
#print(colname)
#function=Group_bet_centrality
#Data=Supra_MST
#Function_output(function,Data,filename,colname,list(range(8)))
#print('we did it')
#--------------------
#MUltilayer_BC_RANDOM

#filename='Bet_cent_FPN_Group_MST_Multi_layer_random_'#+Layer_dic[i]+'_tag_'+str(i)
#colname='Bet_cent_FPN_Group_MST_Multi_layer_random_'#+Layer_dic[i]+'_tag_'+str(i)
#print(filename)
#print(colname)
#function=Group_bet_centrality
#Data=Supra_MST_random
#Function_output(function,Data,filename,colname,list(range(8)))
#print('we did it')
    


#--------------------------------------------------------



#MULTILAYER_ECC
#filename='Eccentricity_FPN_Group_MST_Multi_layer_real_'#+Layer_dic[i]+'_tag_'+str(i)
#colname='Eccentricity_FPN_Group_MST_Multi_layer_real_'#+Layer_dic[i]+'_tag_'+str(i)
#print(filename)
#print(colname)
#function=Group_eccentricity
#Data=Supra_MST
#Function_output(function,Data,filename,colname,list(range(8)))
#print('we did it')

#-----------------
#Multi ecc random
#filename='Eccentricity_FPN_Group_MST_Multi_layer_random_'#+Layer_dic[i]+'_tag_'+str(i)
#colname='Eccentricity_FPN_Group_MST_Multi_layer_random_'#+Layer_dic[i]+'_tag_'+str(i)
#print(filename)
#print(colname)
#function=Group_eccentricity
#Data=Supra_MST_random
#Function_output(function,Data,filename,colname,[i])
#print('we did it')
    
#----------------
#Group_eccentricity_monolayer
    
#for i in range(0,8):
#    filename='Eccentricity_FPN_Group_MST_Mono_layer_real_'+Layer_dic[i]+'_tag_'+str(i)
#    colname='Eccentricity_FPN_Group_MST_Mono_layer_real_'+Layer_dic[i]+'_tag_'+str(i)
#    print(filename)
#    print(colname)
#    function=Group_eccentricity
#    Data=Supra_MST
#    Function_output(function,Data,filename,colname,[i])
#    print('we did it')
#--------------    
# ECC_RANDOM_monolayer

#for i in range(0,8):
#    filename='Eccentricity_FPN_Group_MST_Mono_layer_random_'+Layer_dic[i]+'_tag_'+str(i)
#    colname='Eccentricity_FPN_Group_MST_Mono_layer_random_'+Layer_dic[i]+'_tag_'+str(i)
#    print(filename)
#    print(colname)
#    function=Group_eccentricity
#    Data=Supra_MST_random
#    Function_output(function,Data,filename,colname,[i])
#    print('we did it')
    

# This just get some list - result from a Multilayer computation and convert to the target subnetwork - here is FPN

#Non_norm_Group_eccentricity


#MULTILAYER_ECC - Without Normalization
#filename='Non_norm_Eccentricity_FPN_Group_MST_Multi_layer_real_'#+Layer_dic[i]+'_tag_'+str(i)
#colname='Non_norm_Eccentricity_FPN_Group_MST_Multi_layer_real_'#+Layer_dic[i]+'_tag_'+str(i)
#print(filename)
#print(colname)
#function=Non_norm_Group_eccentricity
#Data=Supra_MST
#Function_output(function,Data,filename,colname,list(range(8)))
#print('we did it')

#-----------------
#Multi ecc random
#filename='Non_norm_Eccentricity_FPN_Group_MST_Multi_layer_random_'#+Layer_dic[i]+'_tag_'+str(i)
#colname='Non_norm_Eccentricity_FPN_Group_MST_Multi_layer_random_'#+Layer_dic[i]+'_tag_'+str(i)
#print(filename)
#print(colname)
#function=Non_norm_Group_eccentricity
#Data=Supra_MST_random
#Function_output(function,Data,filename,colname,[i])
#print('we did it')
    
#----------------
#Group_eccentricity_monolayer
    
#for i in range(0,8):
#    filename='Non_norm_Eccentricity_FPN_Group_MST_Mono_layer_real_'+Layer_dic[i]+'_tag_'+str(i)
#    colname='Non_norm_Eccentricity_FPN_Group_MST_Mono_layer_real_'+Layer_dic[i]+'_tag_'+str(i)
#    print(filename)
#    print(colname)
#    function=Non_norm_Group_eccentricity
#    Data=Supra_MST
#    Function_output(function,Data,filename,colname,[i])
#    print('we did it')
#--------------    
# ECC_RANDOM_monolayer

#for i in range(0,8):
#    filename='Non_norm_Eccentricity_FPN_Group_MST_Mono_layer_random_'+Layer_dic[i]+'_tag_'+str(i)
#    colname='Non_norm_Eccentricity_FPN_Group_MST_Mono_layer_random_'+Layer_dic[i]+'_tag_'+str(i)
#    print(filename)
#    print(colname)
#    function=Non_norm_Group_eccentricity
#    Data=Supra_MST_random
#    Function_output(function,Data,filename,colname,[i])
#    print('we did it')





#MULTILAYER_ECC - Without Normalization and no mask
#filename='Non_norm_Eccentricity_No_Mask_Group_MST_Multi_layer_real_'#+Layer_dic[i]+'_tag_'+str(i)
#colname='Non_norm_Eccentricity_No_Mask_Group_MST_Multi_layer_real_'#+Layer_dic[i]+'_tag_'+str(i)
#print(filename)
#print(colname)
#function=Non_norm_Group_eccentricity
#Data=Supra_MST
#Function_output(function,Data,filename,colname,list(range(8)))
#print('we did it')

#-----------------
#Multi ecc random
#filename='Non_norm_Eccentricity_No_Mask_Group_MST_Multi_layer_random_'#+Layer_dic[i]+'_tag_'+str(i)
#colname='Non_norm_Eccentricity_No_Mask_Group_MST_Multi_layer_random_'#+Layer_dic[i]+'_tag_'+str(i)
#print(filename)
#print(colname)
#function=Non_norm_Group_eccentricity
#Data=Supra_MST_random
#Function_output(function,Data,filename,colname,[i])
#print('we did it')
    
#----------------
#Group_eccentricity_monolayer
    
#for i in range(0,8):
#    filename='Non_norm_Eccentricity_No_Mask_Group_MST_Mono_layer_real_'+Layer_dic[i]+'_tag_'+str(i)
#    colname='Non_norm_Eccentricity_No_Mask_Group_MST_Mono_layer_real_'+Layer_dic[i]+'_tag_'+str(i)
#    print(filename)
#    print(colname)
#    function=Non_norm_Group_eccentricity
#    Data=Supra_MST
#    Function_output(function,Data,filename,colname,[i])
#    print('we did it')
#--------------    
# ECC_RANDOM_monolayer

#for i in range(0,8):
#    filename='Non_norm_Eccentricity_No_Mask_Group_MST_Mono_layer_random_'+Layer_dic[i]+'_tag_'+str(i)
#    colname='Non_norm_Eccentricity_No_Mask_Group_MST_Mono_layer_random_'+Layer_dic[i]+'_tag_'+str(i)
#    print(filename)
#    print(colname)
#   function=Non_norm_Group_eccentricity
#   Data=Supra_MST_random
#    Function_output(function,Data,filename,colname,[i])
#    print('we did it')


#MULTILAYER_ECC - Without Normalization and no mask
#filename='Non_norm_Eccentricity_No_Mask_Group_MST_Multi_layer_real_'#+Layer_dic[i]+'_tag_'+str(i)
#colname='Non_norm_Eccentricity_No_Mask_Group_MST_Multi_layer_real_'#+Layer_dic[i]+'_tag_'+str(i)
#print(filename)
#print(colname)
#function=Non_norm_Group_eccentricity
#Data=Supra_MST
#Function_output(function,Data,filename,colname,list(range(8)))
#print('we did it')

#-----------------
#Multi ecc random
#filename='Non_norm_Eccentricity_No_Mask_Group_MST_Multi_layer_random_'#+Layer_dic[i]+'_tag_'+str(i)
#colname='Non_norm_Eccentricity_No_Mask_Group_MST_Multi_layer_random_'#+Layer_dic[i]+'_tag_'+str(i)
#print(filename)
#print(colname)
#function=Non_norm_Group_eccentricity
#Data=Supra_MST_random
#Function_output(function,Data,filename,colname,[i])
#print('we did it')
    
#----------------
#Group_eccentricity_monolayer
    
#for i in range(0,8):
#    filename='Non_norm_Eccentricity_No_Mask_Group_MST_Mono_layer_real_'+Layer_dic[i]+'_tag_'+str(i)
#    colname='Non_norm_Eccentricity_No_Mask_Group_MST_Mono_layer_real_'+Layer_dic[i]+'_tag_'+str(i)
#    print(filename)
#    print(colname)
#    function=Non_norm_Group_eccentricity
#    Data=Supra_MST
#    Function_output(function,Data,filename,colname,[i])
#    print('we did it')
#--------------    
# ECC_RANDOM_monolayer

#for i in range(0,8):
#    filename='Non_norm_Eccentricity_No_Mask_Group_MST_Mono_layer_random_'+Layer_dic[i]+'_tag_'+str(i)
#    colname='Non_norm_Eccentricity_No_Mask_Group_MST_Mono_layer_random_'+Layer_dic[i]+'_tag_'+str(i)
#    print(filename)
#    print(colname)
#    function=Non_norm_Group_eccentricity
#    Data=Supra_MST_random
#    Function_output(function,Data,filename,colname,[i])
#    print('we did it')


# Same for other stuff
#FPN=[16,17,18,19,20,21,28,29,30,31,93,94,129,130,139,140,169,170]
#Group_bet_cent=Group_bet_centrality(Supra_MST,list(range(13)))
#print('done')
#Group_bet_cent_DWI=Group_bet_centrality(Supra_MST_DWI,list(range(14)))
#Group_bet_cent_FPN=Mask_subnetwork(Group_bet_cent,FPN)
#Group_bet_cent_FPN_DWI=Mask_subnetwork(Group_bet_cent_DWI,FPN)

#SaveSPSS(Group_bet_cent_FPN,'Lucas_Bet_Cen_FPN','BC_Group')
#SaveSPSS(Group_bet_cent_FPN_DWI,'Lucas_Bet_Cen_FPN_DWI','BC_Group_DWI')



#