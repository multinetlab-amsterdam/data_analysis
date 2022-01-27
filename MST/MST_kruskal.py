#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Adaptation of the Kruskal algorithm to generate a maximum spanning tree. 

   The Kruskal algorithm starts with ranking all connection weights from 
   lowest to highest weight. Since we are interested in the strongest
   connections, we ranked all connections from highest to lowest weight
   (formally, our procedure therefore reconstructs the maximum spanning tree).
   We start by disconnecting all nodes, and add the connection with the
   highest weight. Next, the connection with the second highest weight is
   added and this procedure is repeated until all nodes are connected.
   If adding a new connection results in a cycle or loop, this connection is
   discarded, and the next connection ranked by weight is selected.

"""

__author__ = "Mona Zimmermann"
__contact__ = "m.l.m.zimmermann@amsterdamumc.nl"
__date__ = "2021"   ### Date it was created
__status__ = "Concluded" ### Production = still being developed. Else: Concluded/Finished.


####################
# Review History   #
####################

# Reviewed by Lucas Breedt 20210501
# Reviewed by Bernardo Maciel 20220115


####################
# Libraries        #
####################

# Standard imports  #

# Third party imports #
import numpy as np
from tqdm import tqdm


def iscycle(nodes, edge_to_insert):
    """Function to check if setting a link will lead to cycle in graph.

    Reference
    ----------

    Parameters
    ----------
    nodes : 1D numpy array,
        array with length number of nodes present in graph.
        Set of nodes in the graph.

    edge_to_insert : list,
        list containing indices (2) of edge to insert

    Returns
    -------
    nodes: 1D array,
        Array indicating for every node whether it has links or not.
        ´New´ set of nodes.

    c: int,
        indicates if cycle is present or not (0: no cycle, 1: cycle)
    """


    #get highest amount of links + 1
    most_links = max(nodes) + 1

    #c is boolean used to indicate cycle yes or no
    c = False

    #if both nodes dont have links yet
    if (nodes[edge_to_insert[0]] == 0) & (nodes[edge_to_insert[1]] == 0):
        #place links
        nodes[edge_to_insert[0]] = most_links
        nodes[edge_to_insert[1]] = most_links

    #if one node has no links but other does:
    elif nodes[edge_to_insert[0]] == 0:
         #set link to same value as other nodeś link
         nodes[edge_to_insert[0]] = nodes[edge_to_insert[1]]

    #if one node has no links but other does:
    elif nodes[edge_to_insert[1]] == 0:
         #set link to same value as other nodeś link
         nodes[edge_to_insert[1]] = nodes[edge_to_insert[0]]

    #if nodes have same amount of links
    elif nodes[edge_to_insert[0]] == nodes[edge_to_insert[1]]:
        c = True

    #if both nodes have links but not equal amount
    else:
        max_links = max(nodes[edge_to_insert[0]], nodes[edge_to_insert[1]])
        for i in range(0, len(nodes)):
            #if node i has same amount of links as max of edges to insert
            if nodes[i] == max_links:
                #node i is set to lowest amount of links
                nodes[i] = min(nodes[edge_to_insert[0]],nodes[edge_to_insert[1]])

    return(nodes, c)





def kruskal_algorithm(adj_mat):
    """Function to construct the Minnimum Spanning Tree (MST) using the Kruskal algorithm.

    Reference
    ----------


    Parameters
    ----------
    adj_mat : 2D numpy array,
        adjacency matrix containing information on the functional connectivity between regions (size: nr_rois x nr_rois)

    Returns
    -------
    mst: 2D numpy array,
        MST constructed from functional network

    """
    #extract lower triangle of adjacency matrix
    tril_mat = np.tril(adj_mat)

    #find all non-zero elements and place them in new array with shape(3,sth) (1: cols/rows (indx) 2: rows/cols (indx), 3: edge weights)
    PV = np.array((np.where(tril_mat)[0], np.where(tril_mat)[1], tril_mat[np.where(tril_mat)])).T

    #get size of input adjacency matrix --> nr. of nodes
    nr_nodes = adj_mat.shape[0]

    #get size of PV matrix --> nr. of edges (rows)
    nr_edges = PV.shape[0]

    #sort PV by ascending weights order and keep original indices of edges (before sorting)
    sorted_arr = PV[PV[:,2].argsort()]
    
    #resort sorted arr --> strongest weight first
    resorted_arr = np.flipud(sorted_arr)

    #pre-allocate mst matrix
    mst = np.zeros((nr_nodes, nr_nodes))
    nodes = np.zeros(nr_nodes)


    for i in tqdm(range(nr_edges)):
        #edge that will be inserted in graph
        edge_to_insert = resorted_arr[i, [0, 1]]
        edge_to_insert = [int(indx) for indx in edge_to_insert]

        #nodes is array indicating for every node if it has links, c indicates
        #if theres a circle or not when edge is placed
        nodes, c = iscycle(nodes, edge_to_insert)
        if c:
            # if circle created when link is placed --> set rowindex
            # colindex and weight to 0 --> no link will be formed (0 in mst)
            resorted_arr[i,:] = 0

    #create MST adjacency matrix 
    for edge in resorted_arr:
        #if edge is not 0 (row and col indx not 0), then:
        if not (edge[0] == 0 and edge[1] == 0):
            #set edge in resorted_arr to edge weight
            mst[int(edge[0]), int(edge[1])] = edge[2]
            mst[int(edge[1]), int(edge[0])] = edge[2]


    return np.array(mst>0, dtype = 'int'), resorted_arr
