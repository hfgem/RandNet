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
% msgbox('Select folder where the network results are stored.')
% data_path = uigetdir('/Users/hannahgermaine/Documents/PhD/');
% load(strcat(data_path,'/network_spike_sequences_no_I.mat'))

%Store full ranks as matrices
full_ranks = [];
sequence_lengths = [];
for i = 1:length(spike_struct)
    [num_events,~] = size(spike_struct(i).events);
    [n,~] = size(spike_struct(i).spikes_V_m);
    if num_events ~= 0
        for j = 1:num_events
            name = strcat('sequences_',string(j));
            %Update ranks of nonspiking neurons to be num_spiking + 1/2(n - num_spiking)
            full_rank_w_0 = spike_struct(i).ranks.(name);
            ns = length(find(full_rank_w_0));
            new_rank = ns + 0.5*(n - ns);
            new_full_rank = full_rank_w_0;
            new_full_rank(new_full_rank == 0) = new_rank;
            full_ranks = [full_ranks, new_full_rank]; %#ok<AGROW>
            clear full_rank_w_0 ns new_rank new_full_rank
        end  
        clear j num_events
    end
end
clear i   

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