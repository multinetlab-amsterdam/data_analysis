
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 15 11:31:56 2019

@author: turing
"""

"importing important stuff"

""
import numpy as np
import multinetx as mx
import networkx as nx
import matplotlib.pyplot as plt

import numpy as np
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

"This main code is able to create a Multilayer network from a Supra Adjacency Matrix - Those Matrices were made by Lucas - They are both MST and Weighted brain-networks" 


# One should just change the value of (which is now  N=203) as the size of the matrix according to the demand

"You should somehow include the matrices file here - The code is quite robust - It should work with any dataset - If It's prepared well - i. e. - The way Lucas did :-)"

# SOME SETTINGS - Choosing Layer Size and whether you wanna Weighted or Not


#aal_rand_supra_weighted.mat
#aal_supra_weighted_full.mat

layer_size=75 # This are the numbers of nodes in the layer - before was 203
#aal_supra_mst_full.mat

weighted=False

#SHould check which file I should put here!!!!

filename1='aal_supra_mst_full.mat'#'supra_MST_v2'#'supra_MST.mat' #14 layers including dwi
filename2='aal_rand_supra_mst.mat'#'supra_MST_randomized.mat'
#filename3='supra_weighted.mat'
#filename1='supra_MST_full_203rois.mat'
#filename1='supra_MST_first25.mat'
#filename2='supra_first25.mat'
#filename4='supra_weighted_randomized.mat'



# These are the tags for the Multilayer Networks - the way Lucas prepared at the SupraAdjacency Matrix. 
print('0 = fmri, 1=pli delta, 2= pli theta, 3= pli alpha1, 4= pli alpha2, 5 = pli beta, 6 = pli gamma, 7 = DWI.') 

"Loanding the matrices"
#old tags
#print('0 = fmri, 1=pli delta, 2= pli theta, 3= pli alpha1, 4= pli alpha2, 5 = pli beta, 6 = pli gamma, 7 = aec delta, 8 = aec theta,  9 =  aec alpha1, 10 = aec alpha2, 11 = aec beta, 12 = aec gamma, 13 = DWI.') 
#"Loanding the matrices"
#Old_tags
#layer_tags=['0=fmri', '1=pli delta', '2= pli theta', '3= pli alpha1', '4= pli alpha2', '5 = pli beta', '6 = pli gamma', '7 = aec delta','8 = aec theta', '9 =  aec alpha1', '10 = aec alpha2', '11 = aec beta', '12 = aec gamma', '13 = DWI'] 


layer_tags=['0=fmri', '1=pli delta', '2= pli theta', '3= pli alpha1', '4= pli alpha2', '5 = pli beta', '6 = pli gamma', '7 = DWI']
just_tags=['fmri', 'pli_delta', 'pli_theta', 'pli_alpha1', 'pli_alpha2', 'pli_beta', 'pli_gamma', 'DWI'] 
plot_tags=['fMRI', 'PLI delta', 'PLI theta', 'PLI alpha1', 'PLI alpha2', 'PLI beta', 'PLI gamma', 'DWI'] 


#old_dic
#Layer_dic={}
#for i in range(0,14):
#    Layer_dic[i]=just_tags[i]
#print(Layer_dic)
Layer_dic={}
for i in range(0,len(just_tags)):
    Layer_dic[i]=just_tags[i]
print(Layer_dic)


    

#This is the Data in the functions!!! If I can create other data such as random, weighed, etc, I just need to run the special functions!!!!

Supra_MST = scipy.io.loadmat(filename1)
Supra_MST_random =scipy.io.loadmat(filename2)
#Supra_weighted =scipy.io.loadmat(filename3)
#Supra_weighted_random =scipy.io.loadmat(filename4)



#Supra_Weighted=scipy.io.loadmat(filename2)
print(type(Supra_MST)) #Shows that This is in fact a dicionary
#print(type(Supra_Weighted)) #Shows that This is in fact a dicionary

print(Supra_MST.keys())# The command bellow gives the keys
#print(Supra_Weighted.keys())# The command bellow gives the keys

#multilayer=mat['aal_multilayer_adjmats']

#def Create_Multilayer(Data):
    
def Prepare_Multilayer(Data,list_of_layers):
    #"for now should include Supra_MST or Supra_weighted as data and the list of layers you want"
    # I can change for [-1] guy!!
        name=list(Data.keys())[-1]
        multilayer=Data[name]
        
        #multilayer=Data['supra_final']
        where_are_NaNs = np.isnan(multilayer)
        multilayer[where_are_NaNs] = 0
        #number_of_individuals=multilayer.shape[2]
        #layer_numbers=13
        #print(number_of_individuals)
        #layer_numbers=13
        #before it was 205 - make it automatic at some point
        layer_size=75#197 # This are the numbers of nodes in the layer
        N=layer_size
        layer_list=list_of_layers

        layers=[]
        for i in layer_list:
            layers.append(multilayer[(i*N):(i+1)*N,(i*N):(i+1)*N,:])
        return layers


# Just cleaning the data - making nan to zero!!!
# If one includes list(Supra_MST.keys())[-1] might work even not knowing the tag!!!


#This was commented - mayibe need to uncomment

#multilayer=Supra_MST[list(Supra_MST.keys())[-1]]




#multilayer=Supra_MST['supra_final']

#


#where_are_NaNs = np.isnan(multilayer)
#multilayer[where_are_NaNs] = 0


#number_of_individuals=multilayer.shape[2]
#print('The file has '+str(number_of_individuals)+ ' individuals')

"The layers are based in the supra-adjacency matrices - just choose the layers here"

#layer_list=[0,1,2,3,4,5,6,7,8,9,10,11,12]

#
#layer_list=[0,1,2,3,4]




 #Should be 13 here - but it does not converge quick
# THis is the max num of layers




#layer_numbers=13

# Pay attention here could be changed


#layer_size=203 # before was 205
#N=layer_size
#layers=[]
#for i in layer_list:
#    layers.append(multilayer[(i*N):(i+1)*N,(i*N):(i+1)*N,:])


#def I_want_these_layers(layers_you_want):
#    N=203 # before It was 205
#    layers=[]
#    for i in layers_you_want:
#        layers.append(multilayer[(i*N):(i+1)*N,(i*N):(i+1)*N,:])
    #layers
#    return layers
#for i in range(0,layer_numbers):
#    layers.append(multilayer[(i*N):(i+1)*N,(i*N):(i+1)*N,:])
    
# check that, for Linda's data we have!!!
#layer1=multilayer[0:78,0:78,:]
#layer2=multilayer[78:2*78,78:2*78,:]
#layer3=multilayer[2*78:3*78,2*78:3*78,:]



# I Commented until here!!!!


# This function is useful for the agregate - It is the missing step between networkx and Multilayer
    
def agregatte(multiple_layers_list, number_layers):
    k, m = divmod(len(multiple_layers_list), number_layers)
    temp=list(multiple_layers_list[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(number_layers))
    return np.mean(temp, axis=0)

"This creates a multilayer network for each individual"
# This is the new one
def multlayerG(individual,Data,list_of_single_layers):
    "Creates a multilayer for an individual i, knowing the number of layers, and the size of the layers"
    layers= Prepare_Multilayer(Data,list_of_single_layers)
    N =75# 197 # before was 205
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

# These are generic inter connection blocks!!!   
    for i in range(number_of_layers):
        for j in range(number_of_layers):
            if i == j:
                adj_block[i*N:  (i+1)*N,  j*N:(j+1)*N] = np.zeros(N)
            else:
                adj_block[i*N:  (i+1)*N,  j*N:(j+1)*N] = np.identity(N)    # L_12

    mg = mx.MultilayerGraph(list_of_layers=G,inter_adjacency_matrix=adj_block)
    mg.set_edges_weights(intra_layer_edges_weight=1,inter_layer_edges_weight=1)

 
    
    return mg


#def multlayerG(individual,number_of_layers,list_of_layers):
    "Creates a multilayer for an individual i, knowing the number of layers, and the size of the layers"
    

    
# Now we create the group eigenvector centrality for the whole group    


"Here comes all the functions - should do all that we made hypothesis - the way to create the functions is the same"
    
def Group_eingenvector_centrality(Data,list_of_single_layers):

    "This list will save all the eigenvector centralities for all individuals"
    name=list(Data.keys())[-1]

    number_of_individuals=Data[name].shape[2]

    Group_eingenvector=[]
    for individual in range(number_of_individuals):
        temp=multlayerG(individual,Data,list_of_single_layers)

        m=mx.eigenvector_centrality_numpy(temp)
        #m=mx.eigenvector_centrality(multlayerG(individual,number_of_layers,list_of_layers))
    #temp=multlayer3(i)
        temp1=list(m.values())
        temp2=MVagregatte(temp1, len(list_of_single_layers)) 
        #temp2=agregatte(temp1,len(list_of_single_layers))
        # This is a list of lists with all centralities for all individuals
        Group_eingenvector.append(temp2)
        # since we want to buid a flat list 
    flat_list = [item for sublist in Group_eingenvector for item in sublist]
        
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
        temp2=MVagregatte(temp1, len(list_of_single_layers)) 
        #temp2=agregatte(temp1,len(list_of_single_layers))
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
        temp2=temp1#MVagregatte(temp1, len(list_of_single_layers)) 
        
        #temp2=agregatte(temp1,len(list_of_single_layers))
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
        temp2=MVagregatte(temp1, len(list_of_single_layers)) 
        #temp2=agregatte(temp1,len(list_of_single_layers))
        # This is a list of lists with all centralities for all individuals
        Group_bet_centrality.append(temp2)
        # since we want to buid a flat list 
    # Check this flattened = [val for sublist in list_of_lists for val in sublist]
    flat_list = [item for sublist in Group_bet_centrality for item in sublist]
        
    return flat_list



"Here comes the mean and standard deviations from all metrics we analized"
    
def Group_eingenvector_centrality_mean(Data,list_of_single_layers):
    #m=multlayerGNew(3,Supra_MST,[0,5,1])
    "This list will save all the eigenvector centralities means for all individuals"
    name=list(Data.keys())[-1]

    number_of_individuals=Data[name].shape[2]

    Group_eingenvector_mean=[]
    for individual in range(number_of_individuals):
        temp=multlayerG(individual,Data,list_of_single_layers)

        m=mx.eigenvector_centrality_numpy(temp)
        #m=mx.eigenvector_centrality(multlayerG(individual,number_of_layers,list_of_layers))
    #temp=multlayer3(i)
        temp1=list(m.values())
        temp2=MVagregatte(temp1, len(list_of_single_layers)) 
        #temp2=agregatte(temp1,len(list_of_single_layers))
        
        Group_eingenvector_mean.append(np.mean(temp2))
        
    return (Group_eingenvector_mean)





def Group_eingenvector_centrality_std(Data,list_of_single_layers):
    "This list will save all the eigenvector centralities stds for all individuals"
    name=list(Data.keys())[-1]

    number_of_individuals=Data[name].shape[2]
    Group_eingenvector_std=[]
    for individual in range(number_of_individuals):
        temp=multlayerG(individual,Data,list_of_single_layers)
        m=mx.eigenvector_centrality_numpy(temp)
        #m=mx.eigenvector_centrality(multlayerG(individual,number_of_layers,list_of_layers))
    #temp=multlayer3(i)
        temp1=list(m.values())
        temp2=MVagregatte(temp1, len(list_of_single_layers)) # This is MV agregatte - we can change then later for something else
        #temp2=agregatte(temp1,len(list_of_single_layers))
        
        Group_eingenvector_std.append(np.std(temp2))
        
    return (Group_eingenvector_std)



#"This are functions to plot the previous functions"


def Group_hist(Data,list_of_single_layers):
    "this functions aims to plot a histogram with the values of the Eigenvalue centrality for all nodes across all individuals"
    print('layers =',[layer_tags[i] for i in list_of_single_layers])

    temp=Group_eingenvector_centrality(Data,list_of_single_layers)
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

## Include layer tag here! 
## In the end I want to plot all histograms for different twelve layers - the compute the distances and make a matrix

def Plot_eingenvectorcentrality(individual,Data,list_of_single_layers):
    print('layers =',[layer_tags[i] for i in list_of_single_layers])
    #multlayerG(individual,Data,list_of_single_layers)
    m=mx.eigenvector_centrality_numpy(multlayerG(individual,Data,list_of_single_layers))
    #temp=multlayer3(i)
    temp1=list(m.values())
    temp2=MVagregatte(temp1, len(list_of_single_layers)) 
    # This is the Mux Viz Agregatte - We change the agregatte here if yo want later
    #temp2=agregatte(temp1,len(list_of_single_layers))
    plt.hist((temp2))
    # We can edit here the output if we have a vector with the name of the layers
    plt.xlabel('Eigenvector centrality - aggregate ',fontsize=20)
    #plt.xlim(-5,220)
    plt.ylabel("frequence",fontsize=20)
    plt.title('individual '+str(individual))
    plt.show()
    return



def Plot_eingenvectorcentrality_new(individual,Data,list_of_single_layers):
    print('layers =',[layer_tags[i] for i in list_of_single_layers])
    #multlayerG(individual,Data,list_of_single_layers)
    m=mx.eigenvector_centrality_numpy(multlayerG(individual,Data,list_of_single_layers))
    #temp=multlayer3(i)
    temp1=list(m.values())
    temp2=MVagregatte(temp1, len(list_of_single_layers)) 
    # This is the Mux Viz Agregatte - We change the agregatte here if yo want later
    #temp2=agregatte(temp1,len(list_of_single_layers))
    plt.hist((temp2))
    # We can edit here the output if we have a vector with the name of the layers
    plt.xlabel('Eigenvector centrality - aggregate ',fontsize=20)
    #plt.xlim(-5,220)
    plt.ylabel("frequence",fontsize=20)
    plt.ylim(0, 100)
    plt.title('individual '+str(individual))
    plt.show()
    return

def newagregatte(multiple_layers_list, number_layers):
    k, m = divmod(len(multiple_layers_list), number_layers)
    temp=list(multiple_layers_list[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(number_layers))
    #for sublists in temp:
        #m=np.max(temp[sublists])
        #for i in sublists:
            #temp[sublists][i]=temp[sublists][i]/m
            
    return temp#np.mean(temp, axis=0)

# This is the agregatte I suppose Muxviz is using - Double Check with Lucas

def MVagregatte(multiple_layers_list, number_layers):
    k, m = divmod(len(multiple_layers_list), number_layers)
    temp=list(multiple_layers_list[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(number_layers))
    temp_mean=np.mean(temp,axis=0)
    temp_mean=temp_mean/max(temp_mean)
    #for sublists in temp:
    #    m=np.max(temp[sublists])
    #    for i in sublists:
    #        temp[sublists][i]=temp[sublists][i]/m
            
    return temp_mean
            

#
#def agregatte(multiple_layers_list, number_layers):
#    k, m = divmod(len(multiple_layers_list), number_layers)
#    temp=list(multiple_layers_list[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(number_layers))
#    return np.mean(temp, axis=0)

    
#def agregatte(multiple_layers_list, number_layers):
    #k, m = divmod(len(multiple_layers_list), number_layers)
    #temp=list(multiple_layers_list[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(number_layers))
    #return np.mean(temp, axis=0)



def eingenvectorcentrality(individual,Data,list_of_single_layers):
    print('layers =',[layer_tags[i] for i in list_of_single_layers])
    #multlayerG(individual,Data,list_of_single_layers)
    m=mx.eigenvector_centrality_numpy(multlayerG(individual,Data,list_of_single_layers))
    #temp=multlayer3(i)
    temp1=list(m.values())
    #temp2=agregatte(temp1,len(list_of_single_layers))
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
        temp1=list(m.values()) # this is not agregatted 
        temp2=MVagregatte(temp1, len(list_of_single_layers)) # This is Mux Viz Agregatte
        #temp2=agregatte(temp1,len(list_of_single_layers))
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
        temp2=MVagregatte(temp1, len(list_of_single_layers))  # This is Mux Viz agregatte
        #temp2=agregatte(temp1,len(list_of_single_layers)) You can change the agregatte here
        print(temp2)
        Group_degree_centrality_std.append(np.std(temp2))
        
    return (Group_degree_centrality_std)
" We need to create 90 metrics - 45 for the Multilayer and 45 for the random networks"

def Mask_subnetwork(result,target):
    "If you say the target nodes for a given subnetwork, this command you return only the results of list associated with the target nodes"
    chunks = [result[x:x+layer_size] for x in range(0, len(result), layer_size)]
    mask=[chunk[x] for chunk in chunks for x in target]
    
    return mask
# THis creates a SPSS file - Data is the data, name is the name of the file and tag is the collum name
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
    FPN=list(range(0,75)) #THIS IS WITHOUT THE MASK - before 197
    #FPN=[16,17,18,19,20,21,28,29,30,31,93,94,123,124,133,134,163,164]
    #old one - FPN=[16,17,18,19,20,21,28,29,30,31,93,94,129,130,139,140,169,170]
    temp_FPN=Mask_subnetwork(temp,FPN)
    SaveSPSS(temp_FPN,filename,colname)
    
    return

# Create a code - called Running Stuff - So that I can just run things based in this
#Monolayer EC---------------- # Without Mask! Check function above!
#for i in range(0,8):
#    filename='TESTECAAL_No_Mask_Group_MST_Mono_layer_real_'+Layer_dic[i]+'_tag_'+str(i)
#    colname='TESTECAAL_No_Mask_Group_MST_Mono_layer_real_'+Layer_dic[i]+'_tag_'+str(i)
#    print(filename)
#    print(colname)
#    function=Group_eingenvector_centrality
#    Data=Supra_MST
#    Function_output(function,Data,filename,colname,[i])
#    print('we did it')



#M0no RANDOM_EC-----------------
#for i in range(0,8):
#    filename='EC_No_Mask_Group_MST_Mono_layer_random_'+Layer_dic[i]+'_tag_'+str(i)
#    colname='EC_No_Mask_Group_MST_Mono_layer_random_'+Layer_dic[i]+'_tag_'+str(i)
#    print(filename)
#    print(colname)
#    function=Group_eingenvector_centrality
#    Data=Supra_MST_random
#    Function_output(function,Data,filename,colname,[i])
#    print('we did it')
    

#
#MULTILAYER _EC-----real ------
#filename='EC_No_Mask_Group_MST_Multi_layer_real_'#+Layer_dic[i]+'_tag_'+str(i)
#colname='EC_No_Mask_Group_MST_Multi_layer_real_'#+Layer_dic[i]+'_tag_'+str(i)
#print(filename)
#print(colname)
#function=Group_eingenvector_centrality
#Data=Supra_MST
#Function_output(function,Data,filename,colname,list(range(8)))
#print('we did it')

#MULTI_EC------random

#filename='EC_No_Mask_Group_MST_Multi_layer_random_'#+Layer_dic[i]+'_tag_'+str(i)
#colname='EC_No_Mask_Group_MST_Multi_layer_random_'#+Layer_dic[i]+'_tag_'+str(i)
#print(filename)
#print(colname)
#function=Group_eingenvector_centrality
#Data=Supra_MST_random
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