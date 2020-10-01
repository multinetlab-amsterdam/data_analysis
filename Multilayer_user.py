#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  1 14:43:55 2020

@author: fernando
"""

import Multilayer_Main_code as Multinet

# No need to review this file yet- I'm still working in this!


#######################################################################
# HERE ARE NEARLY ALL THE RUNNING ANALYSIS WE DID FOR THE MUMO PROJECT#
#######################################################################


########################
#Monolayer EC - No Mask#
########################

for i in range(0,8):
    filename='Random_No_Mask_Group_MST_Mono_layer_real_'+Layer_dic[i]+'_tag_'+str(i)
    colname='Random_No_Mask_Group_MST_Mono_layer_real_'+Layer_dic[i]+'_tag_'+str(i)
    print(filename)
    print(colname)
    function=Group_eigenvector_centrality
    Data=Supra_MST
    Multinet.Function_output(function,Data,filename,colname,[i])
    print('we did it')


###############
#MULTILAYER-EC#
###############
    
filename='EC_No_Mask_Group_MST_Multi_layer_real_'#+Layer_dic[i]+'_tag_'+str(i)
colname='EC_No_Mask_Group_MST_Multi_layer_real_'#+Layer_dic[i]+'_tag_'+str(i)
print(filename)
print(colname)
function=Group_eigenvector_centrality
Data=Supra_MST
Multinet.Function_output(function,Data,filename,colname,list(range(8)))
print('we did it')


   
##########################
#Monolayer Bet_centrality#
##########################

for i in range(0,8):
    filename='Bet_cent_FPN_Group_MST_Mono_layer_real_'+Layer_dic[i]+'_tag_'+str(i)
    colname='Bet_cent_FPN_Group_MST_Mono_layer_real_'+Layer_dic[i]+'_tag_'+str(i)
    print(filename)
    print(colname)
    function=Group_bet_centrality
    Data=Supra_MST
    Multinet.Function_output(function,Data,filename,colname,[i])
    print('we did it')


################
#MULTILAYER - BC#
################
filename='Bet_cent_FPN_Group_MST_Multi_layer_real_'#+Layer_dic[i]+'_tag_'+str(i)
colname='Bet_cent_FPN_Group_MST_Multi_layer_real_'#+Layer_dic[i]+'_tag_'+str(i)
print(filename)
print(colname)
function=Group_bet_centrality
Data=Supra_MST
Multinet.Function_output(function,Data,filename,colname,list(range(8)))
print('we did it')
#--------------------


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