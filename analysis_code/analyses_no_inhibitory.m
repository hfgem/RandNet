%This script contains code to analyze sequences while excluding inhibitory
%neurons.

%% Load Original Data

%Select the folder where outputs are stored
msgbox('Select folder where the network results are stored.')
data_path = uigetdir('/Users/hannahgermaine/Documents/PhD/');

%Load the membrane potential results, spike events, and parameters
load(strcat(data_path,'/network_var.mat'))
load(strcat(data_path,'/network.mat'))
slashes = find(data_path == '/');
param_path = data_path(1:slashes(end));
load(strcat(param_path,'/parameters.mat'))
clear slashes

%Keep just membrane potential to save space
fields = fieldnames(network_var);
for i = 1:length(fields)
    if ~strcmp(fields{i},'V_m')
        network_var = rmfield(network_var,fields{i});
    end
end
clear i fields


%% Generate Structure Files Excluding Inhibitory Neurons

[network_spike_sequences_no_I, network_cluster_sequences_no_I] = ...
    exclude_inhibitory(network, network_var, parameters);
save(strcat(data_path,'/network_spike_sequences_no_I.mat'),'network_spike_sequences_no_I','-v7.3')
save(strcat(data_path,'/network_cluster_sequences_no_I.mat'),'network_spike_sequences_no_I','-v7.3')

clear network network_var %to save space

%% Calculate Distances if Inhibitory Neurons Were Removed

%Load data if not in workspace
% load(strcat(data_path,'/network_spike_sequences_no_I.mat'))

%Store full ranks as matrices
%Useful Variables
n = parameters.n;

%First reformat ranks into matrix
[full_ranks, sequence_lengths, nonfiring_neurons] = create_rank_matrix(network_spike_sequences_no_I);
num_viable_inits = length(sequence_lengths);
save(strcat(data_path,'/full_ranks_no_I_matrix.mat'),'full_ranks')  

%Calculate distances
full_dist = calculate_vector_distances(full_ranks);
full_dist_vec = nonzeros(triu(full_dist,1));
%Generate shuffled ranks
shuffle_n = 100;
full_shuffle = generate_shuffled_trajectories2(full_ranks, shuffle_n);
%Calculate shuffled distances
full_shuffle_dist = calculate_vector_distances(full_shuffle);
full_shuffle_dist_vec = nonzeros(triu(full_shuffle_dist,1));
%Plot resulting histograms
f = figure;
histogram(full_dist_vec,'DisplayName','Full Vector Distances')
hold on
histogram(full_shuffle_dist_vec,'DisplayName','Shuffled Full Vector Distances')
xlabel('Distance')
ylabel('Number of Distances')
title('Full Rank Sequence Distances, Excluding Inhibitory Neurons')
legend()
savefig(f,strcat(data_path,'/distance_plots_no_I.fig'))
saveas(f,strcat(data_path,'/distance_plots_no_I.jpg'))
saveas(f,strcat(data_path,'/distance_plots_no_I.svg'))