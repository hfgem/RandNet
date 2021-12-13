# RandNet
 This repository contains a collection of code written in MATLAB. The goal of this code is to simulate spike trains using a Leaky Integrate-and-Fire model of firing with a clustered network structure of neurons, with each cluster having randomized sparse connectivity.
 
 ## lif_network_postrotation.m
 This program contains code to initialize parameters, run network creation using create_clusters.m, run a LIF model using lif_sra_calculator_postrotation.m, and perform a number of post-run analyses. The outputs of successful networks are limited by a set of 'success' criteria based on hippocampal research.
 
 ## viualize_cluster_sequences.m
 This program visualizes a sequences of clusters a spike sequences progresses through, based on outputs from lif_network_postrotation.m
 
 ## visualize_voltage_traces.m
 This program visualizes how a spike in one neuron affects connected neurons by looking at membrane potential deflections.
 
 ## STDP_tests.m
 This program runs plasticity tests with a particular network structure and parameter set. The program uses STDP to strengthen one firing initialization, and then has a number of code blocks to visualize and analyze changes in the "learned stimulus" response, as well as other responses.
 
 ## network_tests2.m
 This program runs a grid search of parameter space to determine parameters that result in a network whose behavior meets a strict set of criteria based on biological date. It utilizes the functions parallelize_parameter_tests.m, parallelize_networks.m, and paralellize_network_tests.m to parallelize the grid search and apply the criteria.
 
 ## archive
 This folder contains code that is not currently in use, but was previously useful. It serves as an archive for reference when developing new code.

 ## Functions:
 
 ### calculate_trajectory_similarity_mi.m
 This function calculates the Matching Index (MI) (from Vas et al.) of trajectory similarity between firing sequences. Specifically, it calculates the MI for neuron firing orders including nonspiking neurons at the end of the order, as well as excluding nonspiking neurons.
 
 ### calculate_trajectory_similarity_spearmans.m
 This function calculates the Spearman's Rank Correlation Index of trajectory similarity between firing sequences with unique ranks. Specifically, it calculates the index for neuron firing orders including nonspiking neurons at the end of the order, as well as excluding nonspiking neurons.
 
 ### comp_percentile.m
 This function computes a percentile of a value against a dataset at a 10^(-2) accuracy.
 
 ### create_clusters.m
 This function generates the network clusters and connections based on the number of neurons, number of clusters, number of neurons per cluster, and the probability of connections within a cluster.
 
 ### generate_shuffled_trajectories.m
 This function generates shuffled firing sequences based on the statistics of real trajectories from network simulations.
 
 ### lif_sra_calculator_postrotation.m
 This function uses the leaky integrate-and-fire model of  neuronal firing to calculate the trajectory of membrane potentials, currents, etc... that take place in a particular network with a particular set of parameters and initialization.
 
 ### parallelize_parameter_tests.m
 This function parallelizes tests of different parameter combinations in a grid search and calls parallelize_networks.m to test a number of network initializations for each parameter set. It returns a 'success' value $\in [0,1]$ which represents the fraction of initializations (number of network structures * number of initializations) which successfully pass the response criteria.
 
 ### parallelize_networks.m
 This function parallelizes tests of different network initializations and calls parallelize_network_tests.m to test different firing initializations for each network structure. It returns a 'pass' value $\in [0, num_inits]$ representing the number of initializations that successfully pass the response criteria.
 
 ### paralellize_network_tests.m
 This function initializes firing and then compares the results to a strict set of parameters based on biological literature from Hippocampus studies. It returns a binary 'pass' value of whether or not this initialization passed the response criteria.
 
 ### parallelize_parameter_tests_2.m
 The same as parallelize_parameter_tests.m except it outputs 3 parameters rather than a 'success' value: (1) average number of neurons in a sequence, (2) average firing rate, (3) average event length.
 
 ### parallelize_networks_2.m
 The same as parallelize_networks.m except it outputs 3 parameters rather than a 'success' value: (1) average number of neurons in a sequence, (2) average firing rate, (3) average event length.
 
 ### paralellize_network_tests_2.m
 The same as paralellize_network_tests.m except it outputs 3 parameters rather than a 'success' value: (1) average number of neurons in a sequence, (2) average firing rate, (3) average event length. It also includes a block to generate figures of sequences for successful trials.
 
 
