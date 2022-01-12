# Functions
 This folder contains functions used in the RandNet project.

 ## calculate_trajectory_distances.m
 This function calculates the average distances between mulidimensional trajectories (primarily used for cluster representation trajectories).

 ## calculate_trajectory_similarity_mi.m
 This function calculates the Matching Index (MI) (from Vas et al.) of trajectory similarity between firing sequences. Specifically, it calculates the MI for neuron firing orders including nonspiking neurons at the end of the order, as well as excluding nonspiking neurons.
 
 ## calculate_trajectory_similarity_spearmans.m
 This function calculates the Spearman's Rank Correlation Index of trajectory similarity between firing sequences with unique ranks. Specifically, it calculates the index for neuron firing orders including nonspiking neurons at the end of the order, as well as excluding nonspiking neurons.
 
 ## calculate_vector_distances.m
 This function calculates the distances between vectors representing points in multidimensional space (namely rank vectors).
 
 ## comp_percentile.m
 This function computes a percentile of a value against a dataset at a 10^(-2) accuracy.
 
 ## create_clusters.m
 This function generates the network clusters and connections based on the number of neurons, number of clusters, number of neurons per cluster, and the probability of connections within a cluster.
 
 ## create_rank_matrix.m
 This function pulls ranks from all initializations into a matrix of ranks, a vector of sequence lengths, and a binary matrix of nonspiking neurons.
 
 ## exclude_inhibitory.m
 This function generates a sequence structure file (as in lif_network_postrotation.m) while excluding inhibitory neurons.
 
 ## generate_shuffled_cluster_trajectories.m
This function takes a structure with multi-dimensional trajectories and generates n shuffled trajectories from random sampling.
 
 ## generate_shuffled_trajectories.m
 This function generates shuffled firing sequences based on the statistics of real trajectories from network simulations - uses a structure file.
 
 ## generate_shuffled_trajectories2.m
 This function generates shuffled firing sequences based on the statistics of real trajectories from network simulations - uses a matrix of sequences.
 
 ## lif_sra_calculator_postrotation.m
 This function uses the leaky integrate-and-fire model of  neuronal firing to calculate the trajectory of membrane potentials, currents, etc... that take place in a particular network with a particular set of parameters and initialization.
 
 ## parallelize_parameter_tests_2.m
 The same as parallelize_parameter_tests.m except it outputs 3 parameters rather than a 'success' value: (1) average number of neurons in a sequence, (2) average firing rate, (3) average event length.
 
 ## parallelize_parameter_tests.m
 This function parallelizes tests of different parameter combinations in a grid search and calls parallelize_networks.m to test a number of network initializations for each parameter set. It returns a 'success' value $\in [0,1]$ which represents the fraction of initializations (number of network structures * number of initializations) which successfully pass the response criteria.
 
 ## parallelize_networks_2.m
 The same as parallelize_networks.m except it outputs 3 parameters rather than a 'success' value: (1) average number of neurons in a sequence, (2) average firing rate, (3) average event length.
 
 ## parallelize_networks.m
 This function parallelizes tests of different network initializations and calls parallelize_network_tests.m to test different firing initializations for each network structure. It returns a 'pass' value $\in [0, num_inits]$ representing the number of initializations that successfully pass the response criteria.
 
 ## paralellize_network_tests_2.m
 The same as paralellize_network_tests.m except it outputs 3 parameters rather than a 'success' value: (1) average number of neurons in a sequence, (2) average firing rate, (3) average event length. It also includes a block to generate figures of sequences for successful trials.
 
 ## paralellize_network_tests.m
 This function initializes firing and then compares the results to a strict set of parameters based on biological literature from Hippocampus studies. It returns a binary 'pass' value of whether or not this initialization passed the response criteria.