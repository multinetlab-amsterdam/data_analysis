### Imports ###
import numpy as np
import pandas as pd
from collections import Counter
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

    #get amount of nodes
    nr_nodes = len(nodes)

    #c is boolean used to indicate cycle yes or no
    c = 0

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
        c = 1

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
    sorted_weights = sorted(PV[:,2])#[2].sort()


    old_indx = PV[:,2].argsort()


    #rowindex of sorted weights (so rows(1) = rowindex of weakest weight as sorted in sorted_weights)
    rows = PV[old_indx,0] #[old_indx]


    #colindex of sorted sorted_weights
    cols = PV[old_indx,1] #[old_indx]


    #put in matrix (rows and cols == nodes)
    sorted_arr = np.array((rows, cols, sorted_weights)).T


    #resort sorted arr --> strongest weight first
    resorted_arr = np.flipud(sorted_arr)


    #determine indices of duplicate values --> WHY IS THIS DONE?
    arr = resorted_arr[:,2]

    unique, indx = np.unique(arr, return_index = True)

    duplicate_indx = [value for value in indx if value not in np.arange(len(arr))]


    #set(indx, np.arange(len(arr))) ### CONTINUE

    #pre-allocate mst matrix
    mst = np.zeros((nr_nodes, nr_nodes))
    nodes = np.zeros(nr_nodes)


    for i in tqdm(range(0, nr_edges)):

        #edge that will be inserted in graph
        edge_to_insert = resorted_arr[i, [0, 1]]
        edge_to_insert = [int(indx) for indx in edge_to_insert]



        #nodes is array indicating for every node if it has links, c indicates
        #if theres a circle or not when edge is placed
        nodes, c = iscycle(nodes, edge_to_insert)
        if c == 1:

            #if circle created when link is placed --> set rowindex, colindex and weight to 0 --> no link will be formed (0 in mst)
            resorted_arr[i,:] = 0 #np.zeros(len(resorted_arr[i]))

    #sum of all weights remaining after setting links forming circles to 0
    w = np.sum(resorted_arr[:, 2])

    #create MST adjacency matrix
    for i in range(0, nr_edges):
        #if edge is not 0 (row and col indx not 0), then:
        if (resorted_arr[i, [0,1]] != [0,0]).all():

            #set edge in resorted_arr to edgeweight
            mst[int(resorted_arr[i, 0]), int(resorted_arr[i, 1])] = resorted_arr[i, 2]
            mst[int(resorted_arr[i, 1]), int(resorted_arr[i, 0])] = resorted_arr[i, 2]

        #set edges to larger than 0
        mst = mst>0
        print(mst.shape)
        mst = mst.astype(np.float128)

    return mst
