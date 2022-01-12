# RandNet
 This repository contains a collection of code written in MATLAB. The goal of this code is to simulate spike trains using a Leaky Integrate-and-Fire model of firing with a clustered network structure of neurons, with each cluster having randomized sparse connectivity.
 
 ## lif_network_postrotation.m
 This program contains code to initialize parameters, run network creation using create_clusters.m, run a LIF model using lif_sra_calculator_postrotation.m, and perform a number of post-run analyses. The outputs of successful networks are limited by a set of 'success' criteria based on hippocampal research.
 
 ## network_tests2.m
 This program runs a grid search of parameter space to determine parameters that result in a network whose behavior meets a strict set of criteria based on biological date. It utilizes the functions parallelize_parameter_tests.m, parallelize_networks.m, and paralellize_network_tests.m to parallelize the grid search and apply the criteria.
 
 ## STDP_tests.m
 This program runs plasticity tests with a particular network structure and parameter set. The program uses STDP to strengthen one firing initialization, and then has a number of code blocks to visualize and analyze changes in the "learned stimulus" response, as well as other responses.
 
 ## analysis_code
 This folder contains code dedicated to different analyses of the network outputs.
 
 ## archive
 This folder contains code that is not currently in use, but was previously useful. It serves as an archive for reference when developing new code.
 
 ## functions
 This folder contains functions written for analysis, visualization, reformatting, etc...
