%Calculating the physics of sequences: velocity and momentum

%% Load Data

%Select the folder where outputs are stored
msgbox('Select folder where the network results are stored.')
data_path = uigetdir('/Users/hannahgermaine/Documents/PhD/');

%Load the membrane potential results, spike events, and parameters
load(strcat(data_path,'/network_var.mat'))
load(strcat(data_path,'/network_spike_sequences.mat'))
slashes = find(data_path == '/');
param_path = data_path(1:slashes(end));
load(strcat(param_path,'/parameters.mat'))

%Keep just membrane potential
fields = fieldnames(network_var);
for i = 1:length(fields)
    if ~strcmp(fields{i},'V_m')
        network_var = rmfield(network_var,fields{i});
    end
end
clear i fields

%% Pull Out Data Into New Structure

%Create variable
spike_struct = struct;

%Pull Spikes
for i = 1:length(network_var)
    %Find and store all spike times and neurons
    V_m = network_var(i).V_m;
    spikes_V_m = V_m >= parameters.V_th;
    spike_struct(i).spikes_V_m = spikes_V_m;
    [spikes_x,spikes_t] = find(spikes_V_m);
    spike_times = [spikes_x,spikes_t];
    max_time = max(spikes_t);
    spiking_neurons = unique(spikes_x, 'stable');
    spike_struct(i).spiking_neurons = spiking_neurons;
    spike_struct(i).spike_times = spike_times;
    spike_struct(i).events = network_spike_sequences(i).events;
    %If it's a neuron initialization, store the starting neurons
    if strcmp(parameters.type,'neuron')
        spike_struct(i).init_neur = find(spikes_V_m(:,1));
    end
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

save(strcat(data_path,'/spike_struct.mat'),'spike_struct')
% clear network_var i

%% Analyze New Sequences

%Load data if not in Workspace
% msgbox('Select folder where the network results are stored.')
% data_path = uigetdir('/Users/hannahgermaine/Documents/PhD/');
% load(strcat(data_path,'/spike_struct.mat'))

%Store ranks as matrices
short_ranks = [];
full_ranks = [];
for i = 1:length(spike_struct)
    [num_events,~] = size(spike_struct(i).events);
    if num_events ~= 0
        for j = 1:num_events
            name = strcat('sequences_',string(j));
            full_ranks = [full_ranks, spike_struct(i).ranks.(name)]; %#ok<AGROW>
            short_ranks = [short_ranks, spike_struct(i).short_spike_ranks.(name)]; %#ok<AGROW>
        end  
        clear j num_events
    end
end
clear i    

%Calculate distances
full_dist = calculate_vector_distances(full_ranks);
full_dist_vec = nonzeros(triu(full_dist,1));
short_dist = calculate_vector_distances(short_ranks);
short_dist_vec = nonzeros(triu(short_dist,1));

%Generate shuffled ranks
shuffle_n = 100;
full_shuffle = generate_shuffled_trajectories2(full_ranks, shuffle_n);
short_shuffle = generate_shuffled_trajectories2(short_ranks, shuffle_n);

%Calculate shuffled distances
full_shuffle_dist = calculate_vector_distances(full_shuffle);
full_shuffle_dist_vec = nonzeros(triu(full_shuffle_dist,1));
short_shuffle_dist = calculate_vector_distances(short_shuffle);
short_shuffle_dist_vec = nonzeros(triu(short_shuffle_dist,1));

%Plot resulting histograms
figure;
subplot(1,2,1)
histogram(full_dist_vec,'DisplayName','Full Vector Distances')
hold on
histogram(full_shuffle_dist_vec,'DisplayName','Shuffled Full Vector Distances')
xlabel('Distance')
ylabel('Number of Distances')
title('Full Rank Sequence Distances')
legend()
subplot(1,2,2)
histogram(short_dist_vec,'DisplayName','Short Vector Distances')
hold on
histogram(short_shuffle_dist_vec,'DisplayName','Shuffled Short Vector Distances')
xlabel('Distance')
ylabel('Number of Distances')
title('Short Rank Sequence Distances')
legend()
