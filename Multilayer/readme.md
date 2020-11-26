# Multilayer_Main_code.py - readme
This readme explains how to use the main multilayer code used in the MULTINET lab. This script is used to calculate a multitude of multilayer network metrics: it takes supra-adjacency matrices as input, converts them to multiNetX objects, and defines multiple functions to compute multilayer measures. As most of the code is based on multiNetX, see also the [multiNetX github repository](https://github.com/nkoub/multinetx "multiNetX github repository").

The script consists of three sections. In the first section, the required packages as well as the data are loaded and a quick sanity check is performed. All necessary functions for multilayer analyses are defined in the second section. Finally, in the third section, you will find the main function *function_output* that can be used to calculate any of the multilayer network measures. Note that the majority of the code (i.e. the entirety of Section 2) can thus be executed without any input or modification: if you're only looking to calculate some network measures on  either the provided random data or your own multilayer dataset, please see [Section 1](#section-1:-setting-things-up "Goto Section 1") and [Section 3](#section-3:-calculating-multilayer-metrics") for instructions on how to do so.

##### Table of contents
[Section 1: Setting things up](#section-1-setting-things-up "Goto Section 1")
  * [Input file](#input-file "Goto Input file")  
  * [Required packages](#required-packages "Goto Required packages")  
  * [Settings](#settings "Goto Settings")  
  * [Creating layer tags](#creating-layer-tags "Goto Creating layer tags")  
  * [Loading the matrices](#loading-the-matrices "Goto Loading the matrices")  
  * [Sanity check](#sanity-check "Goto Sanity check")  

[Section 2: Defining the functions](#section-2-defining-the-functions "Goto Section 2")
  * [Preparing the multilayer](#preparing-the-multilayer "Goto Preparing the multilayer")  
  * [Creating the aggregate](#creating-the-aggregate "Goto Creating the aggregate")  
  * [Multilayer functions](#multilayer-functions "Goto Multilayer functions")  
  * [Plotting functions](#plotting-functions "Goto Plotting functions")  
  * [Other functions](#other-functions "Goto Other functions")  

[Section 3: Calculating multilayer metrics](#section-3-calculating-multilayer-metrics")
  * [Function_output](#function_output "Goto Function_output")  

## Section 1: Setting things up
### Input file
As some of the pre-processing of our data was performed using MATLAB and we also constructed our connectivity matrices in MATLAB, this code expects a supra-adjacency matrix of shape regions\*layers *by* regions\*layers *by* subjects, saved in a MATLAB _.mat_ file.

### Required packages
##### Standard imports
* itertools
* multiprocessing _from_ Pool
* time

##### Third party imports
* matplotlib.pyplot _version 3.1.1_
* multinetx
* networkx _version 2.3_
* numpy _version 1.19.0_
* pandas _version 0.25.1_
* pcolor, show, colorbar, xticks, yticks _from_ pylab
* pyreadstat _version 0.2.9_
* scipy _version 1.3.1_
* preprocessing _from_ sklearn _version 0.21.3_
* MinMaxScaler _from_ sklearn.preprocessing

### Settings
Here, the input file and its attributes are defined. 
Specify the number of regions/nodes per layer by changing *layer_size*, and specify whether data is weighted or unweighted by changing the value of *weighted*. Finally, specify the location of the input data under *filename*.

### Creating layer tags
The layer tags created here are used throughout the code to identify the individual layers of the multilayer, and should thus match the ordering of the layers in the supra-adjacency matrix. Change accordingly.

### Loading the matrices
The matlab supra-adjacency matrix specified in [Settings](#settings "Goto Settings") is loaded.

### Sanity check
A quick check to ensure data is loaded correctly. __Note that this is not an automatic process:__ the output should be checked explicitly. Data should be of class 'dict', and the keys of the dictionary should correspond with the layer tags you created in the previous step.

## Section 2: Defining the functions
### Preparing the multilayer
The two functions defined here (*prepare_multilayer* and *multilayer_g*) convert the input data from a Matlab file to a Multinetx friendly object. These functions are used in all the functions for computation of multilayer network measures that are defined later in this code, and do not need any user input.

### Creating the aggregate
This function calculates an aggregate output from nodal multilayer metrics to ensure one value per multilayer node. The aggregate used is the mean of the value of a network property per node in each layer, and is the same method as is used in MuxViz (see the [MuxViz github repository](https://github.com/manlius/muxViz "MuxViz github repository")).

### Multilayer functions
In this section, a multitude of functions for the calculation of multilayer network metrics are defined. A more in-depth description of these functions as well as their required input and output can be found in each functions' docstrings. Briefly, the following functions are defined:
Function | Output
-------- | --------
group_eigenvector_centrality | Nodal eigenvector centrality per subject
group_clustering | Nodal clustering coefficient per subject
group_degree_centrality | Nodal degree centrality per subject
group_eccentricity | Nodal eccentricity per subject
non_norm_group_eccentricity | Non-aggregated nodal eccentricity per subject
group_bet_centrality | Nodal betweenness centrality per subject
group_eigenvector_centrality_mean | Mean eigenvector centrality per subject
group_eigenvector_centrality_std | Standard deviation of eigenvector centrality per subject
eigenvector_centrality | Nodal eigenvector centrality of a single subject
group_degree_centrality_mean | Mean degree centrality per subject
group_degree_centrality_std | Standard deviation of degree centrality per subject

### Plotting functions
Here, two functions to plot histograms of the network metrics are defined:
Function | Output
-------- | --------
plot_group_ec | Histogram of nodal eigenvector centrality across all subjects
plot_ec | Histogram of nodal eigenvector centrality of a single subject

### Other functions
The function *mask_subnetwork* extracts specific nodes from a previously calculated list of network measures (e.g. nodes from the FPN or DMN); *save_csv* saves data (i.e. a list values) to a .csv file for further analysis.

## Section 3: Calculating multilayer metrics
### Function_output
This final function, *function_output*, is the only function that is needed to calculate any of the multilayer network measures described above. It takes five input parameters (described below) and outputs a .csv file containing the desired multilayer network metric.
Input parameter | Description
--------------- | --------------
function        | One of the functions described above
data            | Multilayer data that has been loaded
filename        | Determines name of file to be saved
colname		| Determines column name in savefile
layers		| List of layers to be included for calculation of multilayer measures

