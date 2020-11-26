# Multilayer_Main_code.py - readme
This readme explains how to use the main multilayer code used in the MULTINET lab. This script is used to calculate a multitude of multilayer network metrics: it takes supra-adjacency matrices as input, converts them to multiNetX objects, and defines multiple functions to compute multilayer measures. As most of the code is based on multiNetX, see also the [multiNetX github repository](https://github.com/nkoub/multinetx "multiNetX github repository").

Note that most of the code can be executed without any input or modification. If you're only looking to calculate some network measures on  either the provided random data or your own multilayer dataset, please see the sections [Settings](#Settings "Goto Settings") through [Sanity check](#Sanity check "Goto Sanity check") and [Function_output](#Function_output "Goto Function_output").

##### Table of contents
[Input file](#Input file "Goto Input file")
[Required packages](# Required packages "Goto Required packages")
[Settings](#Settings "Goto Settings")
[Creating layer tags](#Creating layer tags "Goto Creating layer tags")
[Loading the matrices](#Loading the matrices "Goto Loading the matrices")
[Sanity check](#Sanity check "Goto Sanity check")
[Preparing the multilayer](#Preparing the multilayer "Goto Preparing the multilayer")
[Creating the aggregate](#Creating the aggregate "Goto Creating the aggregate")
[Multilayer functions](#Multilayer functions "Goto Multilayer functions")
[Plotting functions](#Plotting functions "Goto Plotting functions")
[Other functions](#Other functions "Goto Other functions")
[Function_output](#Function_output "Goto Function_output")

### Input file
As some of the pre-processing of our data was performed using MATLAB and we also constructed our connectivity matrices in MATLAB, this code expects a supra-adjacency matrix of shape regions\*layers x regions\*layers x subjects, saved in a MATLAB _.mat_ file.

### Required packages
__Standard imports:__
* itertools
* multiprocessing _from_ Pool
* time

__Third party imports:__
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
The matlab supra-adjacency matrix specified in __Settings__ is loaded.

### Sanity check
A quick check to ensure data is loaded correctly. Note that this is not an automatic process: the output should be checked explicitly. Data should be of class 'dict', and the keys of the dictionary should correspond with the layer tags you created in the previous step.

### Preparing the multilayer
The two functions defined here (*Prepare_Multilayer* and *multlayerG*) convert the input data from a Matlab file to a Multinetx friendly object. These functions are used in all the functions for computation of multilayer network measures that are defined later in this code, and do not need any user input.

### Creating the aggregate
This function calculates an aggregate output from nodal multilayer metrics to ensure one value per multilayer node. The aggregate used is the mean of the value of a network property per node in each layer, and is the same method as is used in MuxViz (see [muxviz.net]).

### Multilayer functions
In this section, a multitude of functions for the calculation of multilayer network metrics are defined. A more in-depth description of these functions as well as their required input and output can be found in each functions' docstrings. Briefly, the following functions are defined:
Function | Output
-------- | --------
Group_eigenvector_centrality | Nodal eigenvector centrality per subject
Group_clustering | Nodal clustering coefficient per subject
Group_degree_centrality | Nodal degree centrality per subject
Group_eccentricity | Nodal eccentricity per subject
Non_norm_Group_eccentricity | Non-aggregated nodal eccentricity per subject
Group_bet_centrality | Nodal betweenness centrality per subject
Group_eigenvector_centrality_mean | Mean eigenvector centrality per subject
Group_eigenvector_centrality_std | Standard deviation of eigenvector centrality per subject
Eigenvector_centrality | Nodal eigenvector centrality of a single subject
Group_degree_centrality_mean | Mean degree centrality per subject
Group_degree_centrality_std | Standard deviation of degree centrality per subject

### Plotting functions
Here, two functions to plot histograms of the network metrics are defined:
Function | Output
-------- | --------
Plot_Group_EC | Histogram of nodal eigenvector centrality across all subjects
Plot_EC | Histogram of nodal eigenvector centrality of a single subject

### Other functions
*Mask_subnetwork* can be used to extract specific nodes from a previously calculated list of network measures (e.g. nodes from the FPN or DMN). *SaveSPSS* saves specified data (i.e. a list of values) to a .csv file for further analysis.

#### Function_output
This final function, *Function_output*, is the only function that is needed to calculate any of the multilayer network measures described above. It takes five input parameters (described below) and outputs a .csv file containing the desired multilayer network metric.
Input parameter | Description
--------------- | --------------
function        | One of the functions described above
Data            | Multilayer data that has been loaded
filename        | Determines name of file to be saved
colname		| Determines column name in savefile
layers		| List of layers to be included for calculation of multilayer measures

