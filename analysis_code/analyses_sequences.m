%This script contains code to analyze sequences including inhibitory
%neurons.

saveFlag = 1; % 1 to save analysis results

%% Load Original Data

%Select the folder where outputs are stored
%msgbox('Select folder where the network results are stored.')
data_path = uigetdir('/Users/hannahgermaine/Documents/PhD/');

%Load data
load(strcat(data_path,'/network_spike_sequences.mat'))
load(strcat(data_path,'/network.mat'))
slashes = find(data_path == '/');
param_path = data_path(1:slashes(end));
load(strcat(param_path,'/parameters.mat'))
clear slashes param_path

%Save full_ranks matrix
%First reformat ranks into matrix
[full_ranks, sequence_lengths, nonfiring_neurons] = create_rank_matrix(network_spike_sequences);
num_viable_inits = length(sequence_lengths);
save(strcat(data_path,'/full_ranks_matrix.mat'),'full_ranks')  

%% Basic Stats
%Uses full_ranks matrix and network_spike_sequences
% load(strcat(data_path,'/full_ranks_matrix.mat'))

%First reformat ranks into matrix
[full_ranks, sequence_lengths, nonfiring_neurons] = create_rank_matrix(network_spike_sequences); 
num_viable_inits = length(sequence_lengths);

%Average Length Calculation
avg_length = mean(sequence_lengths);
sd_length = std(sequence_lengths);

%Average Number of Nonfiring Neurons
avg_nonfiring = mean(sum(nonfiring_neurons,1));
frequency_nonfiring = sum(nonfiring_neurons,2)/num_viable_inits;

%Plot nonfiring neuron histogram
figure;
histogram(frequency_nonfiring(frequency_nonfiring ~= 0),10)
xlabel('Frequency of Nonfiring')
ylabel('Number of Neurons')
title('How Frequently Do Nonspiking Neurons Not Spike Across All Sequences?')


%% Sequence Similarity - Spearman's
%Uses network_spike_sequences.mat and parameters.mat
n = parameters.n;
inits = parameters.test_val_max;

%Find which simulated data indices have sequences
viable_inits = [];
for i = 1:inits
    try
        if ~isempty(network_spike_sequences(i).spike_order)
            viable_inits(end+1) = i; %#ok<SAGROW>
        end
    end
end


if length(viable_inits) ~= 1
    %Create a shuffled dataset based on the viable sequences
    shuffle_n = 100;
    %   First generate shuffled data
    [shuffled_spike_sequences] = generate_shuffled_trajectories(n,...
        shuffle_n, network_spike_sequences, viable_inits);
    shuffled_viable_inits = [1:shuffle_n];

    %=====SPEARMAN'S RANK CORRELATIONS=====

    %_____Real Data Ranks_____
    %ranks = Spearman's rank correlation including nonspiking
    %ranks_mod = Spearman's rank correlation excluding nonspiking
    [ranks, ranks_mod] = calculate_trajectory_similarity_spearmans(n, ...
        viable_inits, network_spike_sequences);
    %   Remove NaN Values
    ranks_mod(isnan(ranks_mod)) = 0;
    %   Turn Rank Matrices Into Vectors of Unique Ranks
    ranks_vec = nonzeros(triu(ranks,1)); %all above diagonal
    ranks_mod_vec = nonzeros(triu(ranks_mod,1)); %all above diagonal

    %_____Shuffled Data Ranks_____
    %   Get shuffled ranks
    [shuffled_ranks, shuffled_ranks_mod] = calculate_trajectory_similarity_spearmans(n, ...
        shuffled_viable_inits, shuffled_spike_sequences);
    shuffled_ranks_vec = nonzeros(triu(shuffled_ranks,1)); %all above diagonal
    shuffled_ranks_mod_vec = nonzeros(triu(shuffled_ranks_mod,1)); %all above diagonal

    %_____Calculate Percentiles_____
    ranks_percentile = comp_percentile(shuffled_ranks_vec,mean(ranks_vec));
    ranks_mod_percentile = comp_percentile(shuffled_ranks_mod_vec,mean(ranks_mod_vec));

    %_____Plot Histograms with Percentiles_____
    f = figure;
    ax1 = subplot(1,2,1);
    histogram(shuffled_ranks_vec,'DisplayName','Shuffled Values')
    hold on
    xline(mean(ranks_vec),'-',strcat('Percentile = ',string(ranks_percentile)),'Color','r','DisplayName','Percentile of Average \newline Real Data Value')
    histogram(ranks_vec,'DisplayName','Real Data Values')
    legend()
    title('Spearmans Rank Correlation Rhos with Nonspiking Neurons at End')
    ax2 = subplot(1,2,2);
    histogram(shuffled_ranks_mod_vec,'DisplayName','Shuffled Values')
    hold on
    xline(mean(ranks_mod_vec),'-',strcat('Percentile = ',string(ranks_mod_percentile)),'Color','r','DisplayName','Percentile of Average \newline Real Data Value')
    histogram(ranks_mod_vec,'DisplayName','Real Data Values')
    legend()
    title('Spearmans Rank Correlation Rhos Excluding Nonspiking Neurons')
    f.Position = [187,387,1112,410];
    if saveFlag
        savefig(f,strcat(data_path,'/','spearmans_rank_percentiles.fig'))
        saveas(f,strcat(data_path,'/','spearmans_rank_percentiles.jpg'))
        saveas(f,strcat(data_path,'/','spearmans_rank_percentiles.svg'))
        close(f)
    end
else
    disp('Not Enough Sequences')
end

%% Sequence Similarity - Distance
%Uses network_spike_sequences.mat

%First reformat ranks into matrix
[full_ranks, sequence_lengths, nonfiring_neurons] = create_rank_matrix(network_spike_sequences); 
num_viable_inits = length(sequence_lengths);

%Remove Outliers (uncomment if desired)
%Remove values > 2 sd from mean, < 2 sd from mean, and < 0.05*n
% avg_length = mean(sequence_lengths);
% sd_length = std(sequence_lengths);
% outlier_ind = [find(sequence_lengths > avg_length + 2*sd_length),find(sequence_lengths < avg_length - 2*sd_length)]; %> 2 sd away from mean
% full_ranks_no_outliers = full_ranks;
% full_ranks_no_outliers(:,outlier_ind) = [];

%Calculate distances
full_dist = calculate_vector_distances(full_ranks);
full_dist_vec = nonzeros(triu(full_dist,1));
% full_dist_no_outliers = calculate_vector_distances(full_ranks_no_outliers);
% full_dist_vec_no_outliers = nonzeros(triu(full_dist_no_outliers,1));

%Generate shuffled ranks
shuffle_n = 100;
full_shuffle = generate_shuffled_trajectories2(full_ranks, shuffle_n);
% full_shuffle_no_outliers = generate_shuffled_trajectories2(full_ranks_no_outliers, shuffle_n);

%Calculate shuffled distances
full_shuffle_dist = calculate_vector_distances(full_shuffle);
full_shuffle_dist_vec = nonzeros(triu(full_shuffle_dist,1));
% full_shuffle_dist_no_outliers = calculate_vector_distances(full_shuffle_no_outliers);
% full_shuffle_dist_vec_no_outliers = nonzeros(triu(full_shuffle_dist_no_outliers,1));

%Plot resulting distance histograms
f = figure;
histogram(full_dist_vec,'DisplayName','Full Vector Distances')
hold on
histogram(full_shuffle_dist_vec,'DisplayName','Shuffled Full Vector Distances')
xlabel('Distance')
ylabel('Number of Distances')
title('Full Rank Sequence Distances')
legend()
if saveFlag
    savefig(f,strcat(data_path,'/','dist_hist.fig'))
    saveas(f,strcat(data_path,'/','dist_hist.jpg'))
    saveas(f,strcat(data_path,'/','dist_hist.svg'))
end

% f2 = figure;
% histogram(full_dist_vec_no_outliers,'DisplayName','Full Vector Distances')
% hold on
% histogram(full_shuffle_dist_vec_no_outliers,'DisplayName','Shuffled Full Vector Distances')
% xlabel('Distance')
% ylabel('Number of Distances')
% title({'Full Rank Sequence Distances','No Outliers'})
% legend()
% savefig(f2,strcat(data_path,'/','dist_hist_no_outliers.fig'))
% saveas(f2,strcat(data_path,'/','dist_hist_no_outliers.jpg'))
% saveas(f2,strcat(data_path,'/','dist_hist_no_outliers.svg'))

%% DEPRECATED: Sequence Similarity - Matching Index
%Uses network_spike_sequences.mat and parameters.mat
n = parameters.n;
inits = parameters.test_val_max;

%Find which simulated data indices have sequences
viable_inits = [];
for i = 1:inits
    try
        isempty(network_spike_sequences(i).spike_order);
        if ~isempty(network_spike_sequences(i).spike_order)
            viable_inits(end+1) = i; %#ok<SAGROW>
        end
    end
end

if length(viable_inits) ~= 1
    %Create a shuffled dataset based on the viable sequences
    shuffle_n = 100;
    %   First generate shuffled data
    [shuffled_spike_sequences] = generate_shuffled_trajectories(n,...
        shuffle_n, network_spike_sequences, viable_inits);
    shuffled_viable_inits = [1:shuffle_n];
    %_____Real Data MIs_____
    %ranks = Spearman's rank correlation including nonspiking
    %ranks_mod = Spearman's rank correlation excluding nonspiking
    [matching_index, matching_index_mod] = calculate_trajectory_similarity_mi(n, ...
        viable_inits, network_spike_sequences);
    %   Remove NaN Values
    matching_index_mod(isnan(matching_index_mod)) = 0;
    %   Turn Rank Matrices Into Vectors of Unique Ranks
    matching_index_vec = nonzeros(triu(matching_index,1)); %all above diagonal
    matching_index_mod_vec = nonzeros(triu(matching_index_mod,1)); %all above diagonal

    %_____Shuffled Data MIs_____
    %   Get shuffled ranks
    [shuffled_matching_index, shuffled_matching_index_mod] = calculate_trajectory_similarity_spearmans(n, ...
        shuffled_viable_inits, shuffled_spike_sequences);
    shuffled_matching_index_vec = nonzeros(triu(shuffled_matching_index,1)); %all above diagonal
    shuffled_matching_index_mod_vec = nonzeros(triu(shuffled_matching_index_mod,1)); %all above diagonal

    %_____Calculate Percentiles_____
    matching_index_percentile = comp_percentile(shuffled_matching_index_vec,mean(matching_index_vec));
    matching_index_mod_percentile = comp_percentile(shuffled_matching_index_mod_vec,mean(matching_index_mod_vec));

    %_____Plot Histograms with Percentiles_____
    f = figure;
    ax1 = subplot(1,2,1);
    histogram(shuffled_matching_index_vec,'DisplayName','Shuffled Values')
    hold on
    xline(mean(matching_index_vec),'-',strcat('Percentile = ',string(matching_index_percentile)),'Color','r','DisplayName','Percentile of Average \newline Real Data Value')
    histogram(matching_index_vec,'DisplayName','Real Data Values')
    legend()
    title('Matching Indices with Nonspiking Neurons at End')
    ax2 = subplot(1,2,2);
    histogram(shuffled_matching_index_mod_vec,'DisplayName','Shuffled Values')
    hold on
    try %#ok<TRYNC>
        xline(mean(matching_index_mod_vec),'-',strcat('Percentile = ',string(matching_index_mod_percentile)),'Color','r','DisplayName','Percentile of Average \newline Real Data Value')
        histogram(matching_index_mod_vec,'DisplayName','Real Data Values')
    end
    legend()
    title('Matching Indices Excluding Nonspiking Neurons')
    f.Position = [187,387,1112,410];
    if saveFlag
        savefig(f,strcat(net_save_path,'/','MI_percentiles.fig'))
        saveas(f,strcat(net_save_path,'/','MI_percentiles.jpg'))
        saveas(f,strcat(net_save_path,'/','MI_percentiles.svg'))
        close(f)
    end
else
    disp('Only 1 Sequence')
end
