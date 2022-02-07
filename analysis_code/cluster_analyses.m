%This file contains code blocks dedicated to analyzing the progression of
%spike sequences through different clusters.

%% Load Data

%Select and load specific network data to analyze
net_save_path = uigetdir('/Users/hannahgermaine/Documents/PhD/','Select Network Save Folder'); %Have user input where they'd like the output stored
slashes = find(net_save_path == '/');
save_path = net_save_path(1:slashes(end));
load(strcat(save_path,'/parameters.mat'))
load(strcat(net_save_path,'/network_cluster_sequences.mat'))
load(strcat(net_save_path,'/network_var.mat'))

%Grab relevant information
[~,inits] = size(network_cluster_sequences); %Grab data sizes

%Sequence Types
clust_fields = fieldnames(network_cluster_sequences);
seq_type = clust_fields(3:end);

%Create result save path
cluster_save_path = strcat(net_save_path,'/cluster_analysis/');
if ~isfolder(cluster_save_path)
    mkdir(cluster_save_path);
end

%% Store matrices and structures of cluster ordering

%First we create "rank" vectors for clusters:
%   1. order of first representation (sequence in which a neuron from each
%       cluster first spikes. Note, many neurons are in multiple clusters,
%       so if clusters are represented at the same timepoint, they will be
%       placed in numeral order.
rank_vec = zeros(parameters.clusters,1);
%   2. order of moving sum maximal representation (sequence in which each
%       cluster's most-represented point in the sequence occurs).
rank_sum_vec = zeros(parameters.clusters,1);
%Second we create sequence vectors for clusters:
%   1. the original time-bin sequence is re-written with non-spiking times
%       removed, for use in 1-to-1 trajectory comparison
seq_struct = struct;
%   2. the moving sum sequence is re-written with non-spiking time bins
%       removed, for use in 1-to-1 trajectory comparison
seq_sum_struct = struct;

for s = 1:length(seq_type)
    s_type = seq_type{s};
    if or(strcmp(s_type,'clusters'),strcmp(s_type,'movsum'))
        for i = 1:inits %for each network initialization
            %First check if there are viable events
            if ~isempty(network_cluster_sequences(i).events)
                %Next check that there are cluster sequences
                if ~isempty(network_cluster_sequences(i).clusters)
                    [num_events,~] = size(network_cluster_sequences(i).events);
                    sequences = network_cluster_sequences(i).clusters(1);
                    sequence_names = fieldnames(sequences);
                    for j = 1:num_events %go through all events
                        sequence_j = network_cluster_sequences(i).(s_type).(sequence_names{j});
                        [c_j,t_j] = find(sequence_j);
                        sequence_condensed = sequence_j(:,t_j);
                        if strcmp(s_type,'clusters')
                            j_order = unique(c_j,'stable');
                            ranks_j = zeros(parameters.clusters,1);
                            for k = 1:parameters.clusters
                                k_ind = find(j_order == k);
                                ranks_j(k) = k_ind;
                            end 
                            clear k k_ind
                            sequence_condensed = sequence_j(:,t_j);
                            %Store order of first appearance
                            rank_vec = [rank_vec, ranks_j]; %#ok<*AGROW>
                            %Store cluster sequence condensed
                            seq_struct(end+1).sequence = sequence_condensed; %#ok<*SAGROW>
                        elseif strcmp(s_type,'movsum')
                            max_ind = zeros(parameters.clusters,1);
                            for k = 1:parameters.clusters
                                [~, max_k] = max(sequence_j(k,:));
                                max_ind(k) = max_k;
                            end
                            sort_ind = unique(sort(max_ind));
                            ranks_j = [];
                            for k = 1:length(sort_ind)
                                k_ind = find(max_ind == sort_ind(k));
                                ranks_j = [ranks_j;k_ind];
                            end
                            clear k max_k k_ind max_ind sort_ind
                            %Store order of maximum index first appearance
                            rank_sum_vec = [rank_sum_vec, ranks_j];
                            %Store cluster sequence condensed
                            seq_sum_struct(end+1).sequence = sequence_condensed;
                        end
                        clear c_k t_k j_order ranks_j max_ind 
                    end
                    clear j c_j t_j sequence_j sequence_condensed
                end    
            end
        end
        clear i
    end
end
clear s
%Clear out empty placeholder columns/rows
rank_vec(:,1) = [];
rank_sum_vec(:,1) = [];
seq_struct(1) = [];
seq_sum_struct(1) = [];

%Save results for future use
save(strcat(cluster_save_path,'/rank_vec.mat'),'rank_vec','-v7.3')
save(strcat(cluster_save_path,'/rank_sum_vec.mat'),'rank_sum_vec','-v7.3')
save(strcat(cluster_save_path,'/seq_struct.mat'),'seq_struct','-v7.3')
save(strcat(cluster_save_path,'/seq_sum_struct.mat'),'seq_sum_struct','-v7.3')

%% Calculate cluster sequence similarities

%Load data if not yet in workplace
% load(strcat(cluster_save_path,'/rank_vec.mat'))
% load(strcat(cluster_save_path,'/rank_sum_vec.mat'))
% load(strcat(cluster_save_path,'/seq_struct.mat'))
% load(strcat(cluster_save_path,'/seq_sum_struct.mat'))

%____Description____
%Rank Analyses: we will calculate distances between rank vectors and
%   compare with a shuffled dataset.
%   1. rank_vec = order of initial appearance of a cluster representation,
%       with clusters appearing at the same time in order of index.
%   2. rank_sum_vec = order of maximum firing rate point in sequence.
%Sequence Analyses: we will calculate trajectory similarity by treating
%       each column of the sequence matrix as a point in a 10-dimensional
%       space trajectory, and calculating point-by-point distances and then
%       overall trajectory similarity.
%   1. seq_struct = trajectory similarity of initial timestep data
%   2. seq_sum_struct = trajectory similarity of moving-sum data

%____Rank Analyses____
%Generate shuffles
shuffle_n = 200;
rank_shuffles = generate_shuffled_trajectories2(rank_vec, shuffle_n);
rank_sum_shuffles = generate_shuffled_trajectories2(rank_sum_vec, shuffle_n);
%Calculate true distances
rank_dist = calculate_vector_distances(rank_vec);
rank_sum_dist = calculate_vector_distances(rank_sum_vec);
%Calculate shuffle distances
rank_shuffle_dist = calculate_vector_distances(rank_shuffles);
rank_sum_shuffle_dist = calculate_vector_distances(rank_sum_shuffles);
%Parse out distance vectors
ranks_dist_vec = nonzeros(triu(rank_dist,1));
rank_sum_dist_vec = nonzeros(triu(rank_sum_dist,1));
rank_shuffle_dist_vec = nonzeros(triu(rank_shuffle_dist,1));
rank_sum_shuffle_dist_vec = nonzeros(triu(rank_sum_shuffle_dist,1));
%Plot Histograms of Distributions
f = figure;
ax1 = subplot(1,2,1);
h1 = histogram(ranks_dist_vec);
hold on
h2 = histogram(rank_shuffle_dist_vec);
xlabel('Distance')
ylabel('Number of Values')
title('Rank Distances')
legend('True Rank Distances','Shuffled Rank Distances')
ax2 = subplot(1,2,2);
h3 = histogram(rank_sum_dist_vec);
hold on
h4 = histogram(rank_sum_shuffle_dist_vec);
xlabel('Distance')
ylabel('Number of Values')
title('Maximum Moving Sum Rank Distances')
legend('True Rank Distances','Shuffled Rank Distances')
savefig(f,strcat(cluster_save_path,'/','rank_dist.fig'))
saveas(f,strcat(cluster_save_path,'/','rank_dist.jpg'))
saveas(f,strcat(cluster_save_path,'/','rank_dist.svg'))
clear ax1 h1 h2 ax2 h3 h4

%____Sequence Analyses____
%Generate shuffles
shuffle_n = 200;
seq_shuffles = generate_shuffled_cluster_trajectories(seq_struct, shuffle_n);
seq_sum_shuffles = generate_shuffled_cluster_trajectories(seq_sum_struct, shuffle_n);
%Calculate true distances
[seq_struct_d,seq_struct_d_sd] = calculate_trajectory_distances(seq_struct);
[seq_sum_struct_d,seq_sum_struct_d_sd] = calculate_trajectory_distances(seq_sum_struct);
%Calculate shuffle distances
[seq_shuffle_d,seq_shuffle_d_sd] = calculate_trajectory_distances(seq_shuffles);
[seq_sum_shuffle_d,seq_sum_shuffle_d_sd] = calculate_trajectory_distances(seq_sum_shuffles);
%Parse out distance vectors
seq_struct_d_vec = nonzeros(triu(seq_struct_d,1));
seq_sum_struct_d_vec = nonzeros(triu(seq_sum_struct_d,1));
seq_shuffle_d_vec = nonzeros(triu(seq_shuffle_d,1));
seq_sum_shuffle_d_vec = nonzeros(triu(seq_sum_shuffle_d,1));
%Plot Histograms of Distributions
f2 = figure;
ax1 = subplot(1,2,1);
h1 = histogram(seq_struct_d_vec);
hold on
h2 = histogram(seq_shuffle_d_vec);
xlabel('Distance')
ylabel('Number of Values')
title('Average Trajectory Distances')
legend('True Trajectory Distances','Shuffled Trajectory Distances')
ax2 = subplot(1,2,2);
h3 = histogram(seq_sum_struct_d_vec);
hold on
h4 = histogram(seq_sum_shuffle_d_vec);
xlabel('Distance')
ylabel('Number of Values')
title('Average Moving Sum Trajectory Distances')
legend('True Trajectory Distances','Shuffled Trajectory Distances')
savefig(f2,strcat(cluster_save_path,'/','seq_dist.fig'))
saveas(f2,strcat(cluster_save_path,'/','seq_dist.jpg'))
saveas(f2,strcat(cluster_save_path,'/','seq_dist.svg'))
clear ax1 h1 h2 ax2 h3 h4

%% Are Sequential Neurons In The Same Clusters?
%This block looks at whether successive spikes are primarily from neurons
%in the same clusters or different as compared to shuffled data

% %Select the folder where outputs are stored
% msgbox('Select folder where the network results are stored.')
% data_path = uigetdir('/Users/hannahgermaine/Documents/PhD/');
% 
% %Load Spike Data
% load(strcat(data_path,'/network_spike_sequences_no_I.mat')) %Excluding inhibitory neurons
% load(strcat(data_path,'/network.mat'))

%Store cluster data
cluster_mat = network.cluster_mat;
[n_c,~] = size(cluster_mat); %Number of clusters

%Create Sequence Matrix
[~, ~, ~, sequence_mat] = create_rank_matrix(network_spike_sequences_no_I);
[n,n_s] = size(sequence_mat);

%Determine Sequential Similarity
%This is done by looking at only the neurons which do fire, and their
%sequence, and pairwise calculating how many clusters both neurons belong
%to, and what the average overlap is for the full sequence.
[short_order_sim, short_order] = order_similarity(sequence_mat, cluster_mat);

%Calculate the progression of pair similarity
%Across Trial Averaging
avg_overlap_progression = mean(short_order,2); %Calculate the average overlap at each sequence point
std_overlap_progression = std(short_order,[],2); %Calculate the deviation at each sequence point
curve1 = avg_overlap_progression + std_overlap_progression;
curve2 = avg_overlap_progression - std_overlap_progression;
x = 1:n;
x2 = [x, fliplr(x)];
inBetween = [curve1', fliplr(curve2)'];
%Within Trial Stats
max_block = zeros(1,n_s); %maximum block of overlapping neurons in sequence
for i = 1:n_s %check overlap for each sequence
    zero_s = short_order(:,i) == 0;
    zero_ind = find(zero_s);
    diff = zero_ind(2:end) - zero_ind(1:end-1);
    max_block(1,i) = max(diff);
end 
clear i zero_s zero_ind diff

%Calculate chance similarity
neur_c = sum(cluster_mat,1); %number of clusters each neuron belongs to
avg_neur_c = sum(neur_c)/n; %average number of clusters each neuron belongs to
prob_1_overlap = (avg_neur_c/n_c)^2;
prob_2_overlap = prob_1_overlap*((avg_neur_c - 1)/(n_c-1))^2;

%Null distribution
shuffle_n = 200;
shuffle_seq = generate_shuffled_trajectories2(sequence_mat,shuffle_n);
[short_order_shuffled_sim, short_order_shuffled] = order_similarity(shuffle_seq, cluster_mat);
%Across Trial Averaging
avg_overlap_shuffled_progression = mean(short_order_shuffled,2); %Calculate the average overlap at each sequence point
std_overlap_shuffled_progression = std(short_order_shuffled,[],2); %Calculate the deviation at each sequence point
curve1_shuff = avg_overlap_progression + std_overlap_progression;
curve2_shuff = avg_overlap_progression - std_overlap_progression;
x_shuff = 1:n;
x2_shuff = [x_shuff, fliplr(x_shuff)];
inBetween_shuff = [curve1_shuff', fliplr(curve2_shuff)'];
%Within Trial Stats
max_block_shuffled = zeros(1,n_s); %maximum block of overlapping neurons in sequence
for i = 1:n_s %check overlap for each sequence
    zero_s = short_order_shuffled(:,i) == 0;
    zero_ind = find(zero_s);
    diff = zero_ind(2:end) - zero_ind(1:end-1);
    max_block_shuffled(1,i) = max(diff);
end 
clear i zero_s zero_ind diff

%Plot similarity results
f = figure;
%Histogram of Average Overlap of Each Sequence
subplot(3,1,1)
histogram(short_order_sim,10,'DisplayName','Histogram of Values')
hold on
histogram(short_order_shuffled_sim,10,'DisplayName','Histogram of Null Dist. Values')
xline(prob_1_overlap,'Color','r','LineWidth',1,'DisplayName','Chance Probability of 1 Overlap')
xline(prob_2_overlap,'Color','m','LineWidth',1,'DisplayName','Chance Probability of 2 Overlaps')
xlabel('Pair Cluster Similarity Calculation')
ylabel('Number of Sequences')
title('Histogram of Average Sequence Pair Similarity Values')
legend()
%Histogram of All Overlap Values Across All Sequences
subplot(3,1,2)
histogram(short_order,10,'DisplayName','Histogram of Values')
hold on
histogram(short_order_shuffled,10,'DisplayName','Histogram of Null Dist. Values')
xline(prob_1_overlap,'Color','r','LineWidth',1,'DisplayName','Chance Probability of 1 Overlap')
xline(prob_2_overlap,'Color','m','LineWidth',1,'DisplayName','Chance Probability of 2 Overlaps')
xlabel('Pair Cluster Similarity Calculation')
ylabel('Number of Sequences')
title('Histogram of All Sequence Pair Similarity Values')
legend()
%Average Similarity In Sequence Progression 
subplot(3,1,3)  
%fill(x2, inBetween, 'b','FaceAlpha',0.25,'DisplayName','St. Dev. of Overlap');
hold on
plot(1:n,avg_overlap_progression','LineWidth',1.5,'Color','b','DisplayName','Average Overlap')
%fill(x2_shuff, inBetween_shuff, 'g','FaceAlpha',0.25,'DisplayName','St. Dev. of Null Overlap')
plot(1:n,avg_overlap_shuffled_progression,'LineWidth',1.5,'Color','g','DisplayName','Average Null Overlap')
yline(prob_1_overlap,'Color','r','LineWidth',1,'DisplayName','Chance Probability of 1 Overlap')
yline(prob_2_overlap,'Color','m','LineWidth',1,'DisplayName','Chance Probability of 2 Overlaps')
xlabel('Sequence Index')
ylabel('Cluster Overlap Number')
title('Overlap In Sequence Progression')
legend()
%Save
savefig(f,strcat(data_path,'/cluster_analysis/','sequence_prog_cluster_sim_no_I.fig'))
saveas(f,strcat(data_path,'/cluster_analysis/','sequence_prog_cluster_sim_no_I.jpg'))
saveas(f,strcat(data_path,'/cluster_analysis/','sequence_prog_cluster_sim_no_I.svg'))

%% Percent of Cluster Representation Across Time
%Here we generate matrices representing the percent of representation of
%each cluster at each time point in a sequence, to determine if there is a
%clear progression or trade-off in dominance

% %Select the folder where outputs are stored
% msgbox('Select folder where the network results are stored.')
% data_path = uigetdir('/Users/hannahgermaine/Documents/PhD/');
% 
% %Load Spike Data
% load(strcat(data_path,'/network_spike_sequences_no_I.mat'))
% load(strcat(data_path,'/network.mat'))

%Store Cluster Data
cluster_mat = network.cluster_mat;
cluster_mat_e = cluster_mat(:,network.E_indices);
[n_c,n_e] = size(cluster_mat_e); %Number of clusters

%Calculate percent cluster representation across time
cluster_percent_prog = struct;
for i = 1:length(network_spike_sequences_no_I)
    n_e = network_spike_sequences_no_I(i).num_events;
    if n_e ~= 0
        e_times = network_spike_sequences_no_I(i).events;
        %Cycle through the events
        for j = 1:n_e
            spikes_V_m_e = network_spike_sequences_no_I(i).spikes_V_m(:,e_times(j,1):e_times(j,2));
            clusters_e = cluster_mat_e*spikes_V_m_e;
            %NEEDS MORE WORK
        end    
    end    
    
end   


%% Look at cluster progression from rank vectors

% %Load Spike Data
% load(strcat(data_path,'/network_spike_sequences_no_I.mat'))
% load(strcat(data_path,'/network.mat'))

%Store Cluster Data
cluster_mat = network.cluster_mat;
[n_c,~] = size(cluster_mat); %Number of clusters

%Create Rank Matrix
[full_ranks, sequence_lengths, nonfiring_neurons] = create_rank_matrix(network_spike_sequences);
[n,n_s] = size(full_ranks);

%Store short ranks
short_ranks = zeros(size(full_ranks));
for i = 1:n_s
    short_rank_vec = full_ranks(~nonfiring_neurons(:,i),i);
    short_ranks(1:length(short_rank_vec),i) = short_rank_vec;
end  