%Calculating the physics of sequences: velocity and momentum

%% Load Data

saveFlag = 0; % 1 to save analysis results

%Select the folder where outputs are stored
msgbox('Select folder where the network results are stored.')
data_path = uigetdir('/Users/hannahgermaine/Documents/PhD/');

%Load the membrane potential results, spike events, and parameters
load(strcat(data_path,'/V_m_var.mat'))
load(strcat(data_path,'/network_spike_sequences.mat'))
slashes = find(data_path == '/');
param_path = data_path(1:slashes(end));
load(strcat(param_path,'/parameters.mat'))


%% Pull Out Data Into New Structure

%Create variable
spike_struct = struct;

%Pull Spikes
for i = 1:length(V_m_var)
    %Find and store all spike times and neurons
    V_m = V_m_var(i).V_m;
    spikes_V_m = V_m >= parameters.V_th;
    spike_struct(i).spikes_V_m = spikes_V_m;
    [spikes_x,spikes_t] = find(spikes_V_m);
    spike_times = [spikes_x,spikes_t];
    max_time = max(spikes_t);
    spiking_neurons = unique(spikes_x, 'stable');
    spike_struct(i).spiking_neurons = spiking_neurons;
    spike_struct(i).spike_times = spike_times;
    spike_struct(i).events = network_spike_sequences(i).events;

    %If events are not empty, continue
    if ~isempty(network_spike_sequences(i).events)
        events = network_spike_sequences(i).events;
        [event_num,~] = size(events);
        for j = 1:event_num %Run through each event
            start_ind = find(spike_times(:,2) == events(j,1));
            end_ind = find(spike_times(:,2) == events(j,2));
            event_spike_neurons = spike_times(start_ind:end_ind,1);
            event_spike_times = spike_times(start_ind:end_ind,2);
            neuron_order = unique(event_spike_neurons,'stable');
            ranks = zeros(parameters.n,1);
            for k = 1:parameters.n
                rank_val = find(neuron_order == k,1);
                if ~isempty(rank_val)
                    ranks(k) = rank_val;
                end
            end    
            clear k rank_val
            spike_struct(i).ranks.(strcat('sequences_',string(j))) = ranks;
            spike_struct(i).neuron_order.(strcat('sequences_',string(j))) = neuron_order;
            order_indices = zeros(size(event_spike_neurons)); %y-values
            for k = 1:length(event_spike_neurons)
                ind_k = find(event_spike_neurons(k) == neuron_order);
                order_indices(k) = ind_k;
            end
            clear k ind_k
            order_times = event_spike_times - event_spike_times(1); %x-values
            %Fit a curve to the x and y values to find the point at which
            %the greatest inflection occurs
            [fit_curve, gof] = fit(order_times,order_indices,strcat('-',string(length(neuron_order)),'/(1 + exp(a*(x - b))) + ',string(string(length(neuron_order)))),'StartPoint',[1,max(order_times)/2]);
            [d1,d2] = differentiate(fit_curve,order_times); %first and second derivatives
            [max_d1,max_ind_d1] = max(d1); %Maximum velocity implies greatest inflection point
            [max_d2,max_ind_d2] = max(d2); %Maximum acceleration implies the beginning of the "point of no return" in the sequence
            spike_struct(i).fit_curve.(strcat('sequences_',string(j))) = fit_curve;
            spike_struct(i).gof.(strcat('sequences_',string(j))) = gof;
            spike_struct(i).velocity.(strcat('sequences_',string(j))) = d1;
            spike_struct(i).acceleration.(strcat('sequences_',string(j))) = d2;
            spike_struct(i).max_vel_time.(strcat('sequences_',string(j))) = order_times(max_ind_d1);
            spike_struct(i).max_accel_time.(strcat('sequences_',string(j))) = order_times(max_ind_d2);
            %Store spike sequence up to maximum acceleration point
            spike_struct(i).short_spike_order.(strcat('sequences_',string(j))) = event_spike_neurons(1:max_ind_d2);
            short_ranks = ranks;
            short_ranks(event_spike_neurons(max_ind_d2+1:end)) = 0;
            spike_struct(i).short_spike_ranks.(strcat('sequences_',string(j))) = short_ranks;
            nonspike = ones(parameters.n,1);
            nonspike(event_spike_neurons(1:max_ind_d2)) = 0;
            spike_struct(i).short_nonspiking_neurons.(strcat('sequences_',string(j))) = nonspike;
            clear start_ind end_ind event_spike_neurons event_spike_times ...
                neuron_order ranks order_indices order_times fit_curve d1 ...
                d2 max_d1 max_ind_d1 max_d2 max_ind_d2 nonspike gof short_ranks
        end %end events loop
        clear events event_num j
    end %end empty events loop
    clear V_m spikes_V_m spikes_x spikes_t spike_times max_time spiking_neurons
end   

if saveFlag
    save(strcat(data_path,'/spike_struct.mat'),'spike_struct')
end

%% Analyze Distances in New Sequences

%Load data if not in Workspace
msgbox('Select folder where the network results are stored.')
data_path = uigetdir('/Users/hannahgermaine/Documents/PhD/');
load(strcat(data_path,'/spike_struct.mat'))

%Store ranks as matrices
short_ranks = [];
short_ranks_nonspiking = [];
full_ranks = [];
full_ranks_nonspiking = [];
sequence_lengths = [];
for i = 1:length(spike_struct)
    [num_events,~] = size(spike_struct(i).events);
    [n,~] = size(spike_struct(i).spikes_V_m);
    if num_events ~= 0
        for j = 1:num_events
            name = strcat('sequences_',string(j));
            %Update ranks of nonspiking neurons to be num_spiking + 1/2(n - num_spiking)
            full_rank_w_0 = spike_struct(i).ranks.(name);
            full_ranks_nonspiking = [full_ranks_nonspiking, full_rank_w_0 == 0]; %#ok<AGROW>
            ns = length(find(full_rank_w_0));
            new_rank = ns + 0.5*(n - ns);
            new_full_rank = full_rank_w_0;
            new_full_rank(new_full_rank == 0) = new_rank;
            full_ranks = [full_ranks, new_full_rank]; %#ok<AGROW>
            short_rank_w_0 = spike_struct(i).short_spike_ranks.(name);
            short_ranks_nonspiking = [short_ranks_nonspiking, short_rank_w_0 == 0]; %#ok<AGROW>
            ns = length(find(short_rank_w_0));
            sequence_lengths = [sequence_lengths, ns]; %#ok<AGROW>
            new_rank = ns + 0.5*(n - ns);
            new_short_rank = short_rank_w_0;
            new_short_rank(new_short_rank == 0) = new_rank;
            short_ranks = [short_ranks, new_short_rank]; %#ok<AGROW>
            clear full_rank_w_0 ns new_rank new_full_rank short_rank_w_0 ...
                new_short_rank
        end  
        clear j num_events
    end
end
save(strcat(data_path,'/short_ranks.mat'),'short_ranks')
save(strcat(data_path,'/short_ranks_nonspiking.mat'),'short_ranks_nonspiking')
save(strcat(data_path,'/full_ranks.mat'),'full_ranks')
save(strcat(data_path,'/full_ranks_nonspiking.mat'),'full_ranks_nonspiking')
save(strcat(data_path,'/sequence_lengths.mat'),'sequence_lengths')
clear i   

%Remove outliers
avg_short_length = mean(sequence_lengths);
sd_short_length = std(sequence_lengths);
%Remove values > 2 sd from mean, < 2 sd from mean, and < 0.05*n
outlier_ind = [find(sequence_lengths > avg_short_length + 2*sd_short_length),find(sequence_lengths < avg_short_length - 2*sd_short_length),find(sequence_lengths < 0.05*n)]; %> 2 sd away from mean
full_ranks(:,outlier_ind) = [];
short_ranks(:,outlier_ind) = [];

%Calculate distances
full_dist = calculate_vector_distances(full_ranks);
full_dist_vec = nonzeros(triu(full_dist,1));
short_dist = calculate_vector_distances(short_ranks);
short_dist_vec = nonzeros(triu(short_dist,1));

%Generate shuffled ranks
shuffle_n = 100;
[n,num_seq] = size(full_ranks);
full_shuffle = generate_shuffled_trajectories2(full_ranks, shuffle_n);
short_shuffle = generate_shuffled_trajectories2(short_ranks, shuffle_n);
short_shuffle_lengths = [];
for i = 1:shuffle_n
    [mode_i, mode_i_freq] = mode(short_shuffle(:,i));
    short_shuffle_lengths = [short_shuffle_lengths,n - mode_i_freq];
end    

%Calculate shuffled distances
full_shuffle_dist = calculate_vector_distances(full_shuffle);
full_shuffle_dist_vec = nonzeros(triu(full_shuffle_dist,1));
short_shuffle_dist = calculate_vector_distances(short_shuffle);
short_shuffle_dist_vec = nonzeros(triu(short_shuffle_dist,1));

%Plot resulting histograms
figure;
subplot(1,2,1)
histogram(full_shuffle_dist_vec,'DisplayName','Shuffled Full Vector Distances')
hold on
histogram(full_dist_vec,'DisplayName','Full Vector Distances')
xlabel('Distance')
ylabel('Number of Distances')
title('Full Rank Sequence Distances')
legend()
subplot(1,2,2)
histogram(short_shuffle_dist_vec,'DisplayName','Shuffled Short Vector Distances')
hold on
histogram(short_dist_vec,'DisplayName','Short Vector Distances')
xlabel('Distance')
ylabel('Number of Distances')
title('Short Rank Sequence Distances')
legend()

%% Matching Indices with Fractional Rank Nonspiking

%Load data if not in Workspace
msgbox('Select folder where the network results are stored.')
data_path = uigetdir('/Users/hannahgermaine/Documents/PhD/');
load(strcat(data_path,'/spike_struct.mat')) %From Above Code Block
load(strcat(data_path,'/short_ranks.mat')) %From Above Code Block
load(strcat(data_path,'/short_ranks_nonspiking.mat')) %From Above Code Block
load(strcat(data_path,'/full_ranks.mat')) %From Above Code Block
load(strcat(data_path,'/full_ranks_nonspiking.mat')) %From Above Code Block
load(strcat(data_path,'/sequence_lengths.mat')) %From Above Code Block

%___Calculate Matching Indices (Vas et al..)___
shuffle_n = 100;

%Generate shuffled sequences with fractional ranks
full_shuffle_spike = generate_shuffled_trajectories2(full_ranks, shuffle_n);
full_shuffle_nonspike = [];
for i = 1:shuffle_n
    full_shuffle_spike_seq = full_shuffle_spike(:,i);
    [mode_full_shuffle_spike_seq, mode_n] = mode(full_shuffle_spike_seq);
    if mode_n > 1
        full_shuffle_nonspike = [full_shuffle_nonspike, full_shuffle_spike_seq == mode_full_shuffle_spike_seq];
    else
        full_shuffle_nonspike = [full_shuffle_nonspike, zeros(size(full_shuffle_spike_seq))];
    end
end    
short_shuffle_spike = generate_shuffled_trajectories2(short_ranks, shuffle_n);
short_shuffle_nonspike = [];
for i = 1:shuffle_n
    short_shuffle_spike_seq = short_shuffle_spike(:,i);
    [mode_short_shuffle_spike_seq, mode_n] = mode(short_shuffle_spike_seq);
    if mode_n > 1
        short_shuffle_nonspike = [short_shuffle_nonspike, short_shuffle_spike_seq == mode_short_shuffle_spike_seq];
    else
        short_shuffle_nonspike = [short_shuffle_nonspike, zeros(size(short_shuffle_spike_seq))];
    end
end
%Matching indices with ranks of nonspiking = n_s + (n - n_s)/2
penalize_nonspike = 1;
[full_matching_index, full_matching_index_mod] = calculate_trajectory_similarity_mi2(full_ranks, ...
    full_ranks_nonspiking, penalize_nonspike);
[short_matching_index, short_matching_index_mod] = calculate_trajectory_similarity_mi2(short_ranks, ...
    short_ranks_nonspiking, penalize_nonspike);
[shuffle_full_matching_index, shuffle_full_matching_index_mod] = calculate_trajectory_similarity_mi2(full_shuffle_spike, ...
    full_shuffle_nonspike, penalize_nonspike);
[shuffle_short_matching_index, shuffle_short_matching_index_mod] = calculate_trajectory_similarity_mi2(short_ranks, ...
    short_shuffle_nonspike, penalize_nonspike);
%Vectorize the MI Values
full_matching_index_vec = nonzeros(triu(full_matching_index,1));
full_matching_index_mod_vec = nonzeros(triu(full_matching_index_mod,1));
short_matching_index_vec = nonzeros(triu(short_matching_index,1));
short_matching_index_mod_vec = nonzeros(triu(short_matching_index_mod,1));
shuffle_full_matching_index_vec = nonzeros(triu(shuffle_full_matching_index,1));
shuffle_full_matching_index_mod_vec = nonzeros(triu(shuffle_full_matching_index_mod,1));
shuffle_short_matching_index_vec = nonzeros(triu(shuffle_short_matching_index,1));
shuffle_short_matching_index_mod_vec = nonzeros(triu(shuffle_short_matching_index,1));
%Plot Histograms of MI Values for fractional ranks
figure;
subplot(2,2,1)
histogram(shuffle_full_matching_index_vec,'DisplayName','Shuffled Full Sequence Matching Index')
hold on
histogram(full_matching_index_vec,'DisplayName','Full Sequence Matching Index')
xlabel('Matching Index')
ylabel('Number of Occurrences')
title('Full Sequence Matching Indices')
legend()
subplot(2,2,2)
histogram(shuffle_short_matching_index_vec,'DisplayName','Shuffled Short Sequence Matching Index')
hold on
histogram(short_matching_index_vec,'DisplayName','Short Sequence Matching Index')
xlabel('Matching Index')
ylabel('Number of Occurrences')
title('Short Sequence Matching Indices')
legend()
subplot(2,2,3)
histogram(shuffle_full_matching_index_mod_vec,'DisplayName','Shuffled Full Sequence Matching Index')
hold on
histogram(full_matching_index_mod_vec,'DisplayName','Full Sequence Matching Index')
xlabel('Matching Index')
ylabel('Number of Occurrences')
title({'Full Sequence Matching Indices','Excluding Overlapping Nonspiking'})
legend()
subplot(2,2,4)
histogram(shuffle_short_matching_index_mod_vec,'DisplayName','Shuffled Short Sequence Matching Index')
hold on
histogram(short_matching_index_mod_vec,'DisplayName','Short Sequence Matching Index')
xlabel('Matching Index')
ylabel('Number of Occurrences')
title({'Short Sequence Matching Indices','Excluding Overlapping Nonspiking'})
legend()
%Z-Scores of MI Values
full_MI_z_score = (full_matching_index_vec - mean(full_matching_index_vec))/std(full_matching_index_vec);
full_MI_mod_z_score = (full_matching_index_mod_vec - mean(full_matching_index_mod_vec))/std(full_matching_index_mod_vec);
short_MI_z_score = (short_matching_index_vec - mean(short_matching_index_vec))/std(short_matching_index_vec);
short_MI_mod_z_score = (short_matching_index_mod_vec - mean(short_matching_index_mod_vec))/std(short_matching_index_mod_vec);
shuffle_full_MI_z_score = (shuffle_full_matching_index_vec - mean(shuffle_full_matching_index_vec))/std(shuffle_full_matching_index_vec);
shuffle_full_MI_mod_z_score = (shuffle_full_matching_index_mod_vec - mean(shuffle_full_matching_index_mod_vec))/std(shuffle_full_matching_index_mod_vec);
shuffle_short_MI_z_score = (shuffle_short_matching_index_vec - mean(shuffle_short_matching_index_vec))/std(shuffle_short_matching_index_vec);
shuffle_short_MI_mod_z_score = (shuffle_short_matching_index_mod_vec - mean(shuffle_short_matching_index_mod_vec))/std(shuffle_short_matching_index_mod_vec); 
%Plot Histograms of Z Scores
figure;
subplot(2,2,1)
histogram(shuffle_full_MI_z_score,'DisplayName','Shuffled')
hold on
histogram(full_MI_z_score,'DisplayName','Full Sequence')
xlabel('Matching Index')
ylabel('Number of Occurrences')
title('Full Sequence Matching Index Z-Scores')
legend()
subplot(2,2,2)
histogram(shuffle_short_MI_z_score,'DisplayName','Shuffled')
hold on
histogram(short_MI_z_score,'DisplayName','Short Sequence')
xlabel('Matching Index')
ylabel('Number of Occurrences')
title('Short Sequence Matching Index Z-Scores')
legend()
subplot(2,2,3)
histogram(shuffle_full_MI_mod_z_score,'DisplayName','Shuffled')
hold on
histogram(full_MI_mod_z_score,'DisplayName','Full Sequence')
xlabel('Matching Index')
ylabel('Number of Occurrences')
title({'Full Sequence Matching Index Z-Scores','Excluding Overlapping Nonspiking'})
legend()
subplot(2,2,4)
histogram(shuffle_short_MI_mod_z_score,'DisplayName','Shuffled')
hold on
histogram(short_MI_mod_z_score,'DisplayName','Short Sequence')
xlabel('Matching Index')
ylabel('Number of Occurrences')
title({'Short Sequence Matching Index Z-Scores','Excluding Overlapping Nonspiking'})
legend()

%% Matching Index with 0 Rank Nonspiking

%Uncomment next section if running independently of the above MI fractional
%rank section.

% %Load data if not in Workspace
% msgbox('Select folder where the network results are stored.')
% data_path = uigetdir('/Users/hannahgermaine/Documents/PhD/');
% load(strcat(data_path,'/spike_struct.mat')) %From Above Code Block
% load(strcat(data_path,'/short_ranks.mat')) %From Above Code Block
% load(strcat(data_path,'/short_ranks_nonspiking.mat')) %From Above Code Block
% load(strcat(data_path,'/full_ranks.mat')) %From Above Code Block
% load(strcat(data_path,'/full_ranks_nonspiking.mat')) %From Above Code Block
% load(strcat(data_path,'/sequence_lengths.mat')) %From Above Code Block
% 
% %___Calculate Matching Indices (Vas et al..)___
% shuffle_n = 100;
% 
% %Generate shuffled sequences with fractional ranks
% full_shuffle_spike = generate_shuffled_trajectories2(full_ranks, shuffle_n);
% full_shuffle_nonspike = [];
% for i = 1:shuffle_n
%     full_shuffle_spike_seq = full_shuffle_spike(:,i);
%     [mode_full_shuffle_spike_seq, mode_n] = mode(full_shuffle_spike_seq);
%     if mode_n > 1
%         full_shuffle_nonspike = [full_shuffle_nonspike, full_shuffle_spike_seq == mode_full_shuffle_spike_seq];
%     else
%         full_shuffle_nonspike = [full_shuffle_nonspike, zeros(size(full_shuffle_spike_seq))];
%     end
% end    
% short_shuffle_spike = generate_shuffled_trajectories2(short_ranks, shuffle_n);
% short_shuffle_nonspike = [];
% for i = 1:shuffle_n
%     short_shuffle_spike_seq = short_shuffle_spike(:,i);
%     [mode_short_shuffle_spike_seq, mode_n] = mode(short_shuffle_spike_seq);
%     if mode_n > 1
%         short_shuffle_nonspike = [short_shuffle_nonspike, short_shuffle_spike_seq == mode_short_shuffle_spike_seq];
%     else
%         short_shuffle_nonspike = [short_shuffle_nonspike, zeros(size(short_shuffle_spike_seq))];
%     end
% end

%Generate shuffled (and true) sequences with 0 ranks
full_ranks_0_nonspike = full_ranks.*(~full_ranks_nonspiking);
short_ranks_0_nonspike = short_ranks.*(~short_ranks_nonspiking);
full_shuffle_0_nonspike = full_shuffle_spike.*(~full_shuffle_nonspike);
short_shuffle_0_nonspike = short_shuffle_spike.*(~short_shuffle_nonspike);
%Matching indices with ranks of nonspiking = 0
penalize_nonspike = 1;
[full_0_matching_index, full_0_matching_index_mod] = calculate_trajectory_similarity_mi2(full_ranks_0_nonspike, ...
    full_ranks_nonspiking, penalize_nonspike);
[short_0_matching_index, short_0_matching_index_mod] = calculate_trajectory_similarity_mi2(short_ranks_0_nonspike, ...
    short_ranks_nonspiking, penalize_nonspike);
[shuffle_full_0_matching_index, shuffle_full_0_matching_index_mod] = calculate_trajectory_similarity_mi2(full_shuffle_0_nonspike, ...
    full_shuffle_nonspike, penalize_nonspike);
[shuffle_short_0_matching_index, shuffle_short_0_matching_index_mod] = calculate_trajectory_similarity_mi2(short_shuffle_0_nonspike, ...
    short_shuffle_nonspike, penalize_nonspike);
%Vectorized MI Values
full_0_matching_index_vec = nonzeros(triu(full_0_matching_index,1));
full_0_matching_index_mod_vec = nonzeros(triu(full_0_matching_index_mod,1));
short_0_matching_index_vec = nonzeros(triu(short_0_matching_index,1));
short_0_matching_index_mod_vec = nonzeros(triu(short_0_matching_index_mod,1));
shuffle_full_0_matching_index_vec = nonzeros(triu(shuffle_full_0_matching_index,1));
shuffle_full_0_matching_index_mod_vec = nonzeros(triu(shuffle_full_0_matching_index_mod,1));
shuffle_short_0_matching_index_vec = nonzeros(triu(shuffle_short_0_matching_index,1));
shuffle_short_0_matching_index_mod_vec = nonzeros(triu(shuffle_short_0_matching_index_mod,1));
%Plot Histograms of MI Values
figure;
subplot(2,2,1)
histogram(shuffle_full_0_matching_index_vec,'DisplayName','Shuffled Full Sequence Matching Index')
hold on
histogram(full_0_matching_index_vec,'DisplayName','Full Sequence Matching Index')
xlabel('Matching Index')
ylabel('Number of Occurrences')
title('Full Sequence Matching Indices')
legend()
subplot(2,2,2)
histogram(shuffle_short_0_matching_index_vec,'DisplayName','Shuffled Short Sequence Matching Index')
hold on
histogram(short_0_matching_index_vec,'DisplayName','Short Sequence Matching Index')
xlabel('Matching Index')
ylabel('Number of Occurrences')
title('Short Sequence Matching Indices')
legend()
subplot(2,2,3)
histogram(shuffle_full_0_matching_index_mod_vec,'DisplayName','Shuffled Full Sequence Matching Index')
hold on
histogram(full_0_matching_index_mod_vec,'DisplayName','Full Sequence Matching Index')
xlabel('Matching Index')
ylabel('Number of Occurrences')
title({'Full Sequence Matching Indices','Excluding Overlapping Nonspiking'})
legend()
subplot(2,2,4)
histogram(shuffle_short_0_matching_index_mod_vec,'DisplayName','Shuffled Short Sequence Matching Index')
hold on
histogram(short_0_matching_index_mod_vec,'DisplayName','Short Sequence Matching Index')
xlabel('Matching Index')
ylabel('Number of Occurrences')
title({'Short Sequence Matching Indices','Excluding Overlapping Nonspiking'})
legend()
%Z-Scores of MI Values
full_0_MI_z_score = (full_0_matching_index_vec - mean(full_0_matching_index_vec))/std(full_0_matching_index_vec);
full_0_MI_mod_z_score = (full_0_matching_index_mod_vec - mean(full_0_matching_index_mod_vec))/std(full_0_matching_index_mod_vec);
short_0_MI_z_score = (short_0_matching_index_vec - mean(short_0_matching_index_vec))/std(short_0_matching_index_vec);
short_0_MI_mod_z_score = (short_0_matching_index_mod_vec - mean(short_0_matching_index_mod_vec))/std(short_0_matching_index_mod_vec);
shuffle_full_0_MI_z_score = (shuffle_full_0_matching_index_vec - mean(shuffle_full_0_matching_index_vec))/std(shuffle_full_0_matching_index_vec);
shuffle_full_0_MI_mod_z_score = (shuffle_full_0_matching_index_mod_vec - mean(shuffle_full_0_matching_index_mod_vec))/std(shuffle_full_0_matching_index_mod_vec);
shuffle_short_0_MI_z_score = (shuffle_short_0_matching_index_vec - mean(shuffle_short_0_matching_index_vec))/std(shuffle_short_0_matching_index_vec);
shuffle_short_0_MI_mod_z_score = (shuffle_short_0_matching_index_mod_vec - mean(shuffle_short_0_matching_index_mod_vec))/std(shuffle_short_0_matching_index_mod_vec); 
%Plot Histograms of Z Scores
figure;
subplot(2,2,1)
histogram(shuffle_full_0_MI_z_score,'DisplayName','Shuffled')
hold on
histogram(full_0_MI_z_score,'DisplayName','Full Sequence')
xlabel('Matching Index')
ylabel('Number of Occurrences')
title('Full Sequence Matching Index Z-Scores')
legend()
subplot(2,2,2)
histogram(shuffle_short_0_MI_z_score,'DisplayName','Shuffled')
hold on
histogram(short_0_MI_z_score,'DisplayName','Short Sequence')
xlabel('Matching Index')
ylabel('Number of Occurrences')
title('Short Sequence Matching Index Z-Scores')
legend()
subplot(2,2,3)
histogram(shuffle_full_0_MI_mod_z_score,'DisplayName','Shuffled')
hold on
histogram(full_0_MI_mod_z_score,'DisplayName','Full Sequence')
xlabel('Matching Index')
ylabel('Number of Occurrences')
title({'Full Sequence Matching Index Z-Scores','Excluding Overlapping Nonspiking'})
legend()
subplot(2,2,4)
histogram(shuffle_short_0_MI_mod_z_score,'DisplayName','Shuffled')
hold on
histogram(short_0_MI_mod_z_score,'DisplayName','Short Sequence')
xlabel('Matching Index')
ylabel('Number of Occurrences')
title({'Short Sequence Matching Index Z-Scores','Excluding Overlapping Nonspiking'})
legend()

%% Investigating interesting results
%Full Rank current sequences exhibited a bimodal distance distribution -
%investigating the differences in the peaks by splitting at a distance of
%500.

[n,~] = size(full_ranks);
[row, col] = find(full_dist < 200);
close_pairs = [0,0];
close_pair_lengths = [0,0];
for i = 1:length(row)
    i_ind = row(i);
    j_ind = col(i);
    if i_ind ~= j_ind
        exists_already = find(sum(close_pairs - ones(size(close_pairs)).*[j_ind,i_ind],2) == 0,1);
        if isempty(exists_already)
            close_pairs = [close_pairs; i_ind, j_ind];
            [mode_i, mode_i_freq] = mode(full_ranks(:,i_ind));
            i_length = n - mode_i_freq;
            [mode_j, mode_j_freq] = mode(full_ranks(:,j_ind));
            j_length = n - mode_j_freq;
            close_pair_lengths = [close_pair_lengths; i_length, j_length];
        end
    end
end    
close_pairs(1,:) = [];
close_pair_lengths(1,:) = [];
close_length_diff = abs(close_pair_lengths(:,2) - close_pair_lengths(:,1));
close_unique = unique(close_pairs);
close_unique_lengths = sum(mod(full_ranks(:,close_unique),1) == 0,1);
avg_close_unique_length = mean(close_unique_lengths);


[row, col] = find(full_dist >= 500);
far_pairs = [0,0];
far_pair_lengths = [0,0];
for i = 1:length(row)
    i_ind = row(i);
    j_ind = col(i);
    if i_ind ~= j_ind
        exists_already = find(sum(far_pairs - ones(size(far_pairs)).*[j_ind,i_ind],2) == 0,1);
        if isempty(exists_already)
            far_pairs = [far_pairs; i_ind, j_ind];
            [mode_i, mode_i_freq] = mode(full_ranks(:,i_ind));
            i_length = n - mode_i_freq;
            [mode_j, mode_j_freq] = mode(full_ranks(:,j_ind));
            j_length = n - mode_j_freq;
            far_pair_lengths = [far_pair_lengths; i_length, j_length];
        end
    end
end
far_pairs(1,:) = [];
far_pair_lengths(1,:) = [];
far_length_diff = abs(far_pair_lengths(:,2) - far_pair_lengths(:,1));
far_unique = unique(far_pairs);
far_unique_lengths = sum(mod(full_ranks(:,far_unique),1) == 0,1);
avg_far_unique_length = mean(far_unique_lengths);

figure;
h1 = histogram(close_length_diff);
hold on
h2 = histogram(far_length_diff);
xlabel('Distance Difference')
ylabel('Number of Sequence Pairs')
legend('Close Sequence Length Differences', 'Far Sequence Length Differences')

%Length distribution ks test
[h_len,p_len] = kstest2(close_unique_lengths,far_unique_lengths);
%Difference in length ks test
[h_diff,p_diff] = kstest2(close_length_diff,far_length_diff);


