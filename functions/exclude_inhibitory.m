%Generate sequence file excluding inhibitory neurons.

function [network_spike_sequences_no_I, network_cluster_sequences_no_I] = ...
    exclude_inhibitory(network, network_var, parameters)
    %ABOUT: This function takes structure outputs from
    %lif_network_postrotation.m and generates new structures with the same
    %data except calculated with inhibitory neurons excluded.
    %
    %INPUTS:
    %   1. network = a structure containing high level network structure
    %       information, including which neurons belong to which clusters.
    %   2. network_var = a structure containing the membrane potential
    %       resulting from the original simulation.
    %   3. parameters = a structure containing parameters used in the
    %       original simulation.
    %OUTPUTS:
    %   1. network_spike_sequences_no_I = a structure containing events and
    %       their corresponding rank and sequence data.
    %   2. network_cluster_sequences_no_I = a structure containing events
    %       and their corresponding cluster sequence data.
    
    e_indices = network.E_indices;
    %Set up structures
    network_spike_sequences_no_I = struct;
    network_cluster_sequences_no_I = struct;
    %Pull Spikes
    for i = 1:length(network_var)
        %Find and store all spike times, excluding inhibitory neurons
        V_m = network_var(i).V_m;
        spikes_V_m = V_m(e_indices,:) >= parameters.V_th;
        network_spike_sequences_no_I(i).spikes_V_m = spikes_V_m;
        %Find events
        [spikes_x,spikes_t] = find(spikes_V_m);
        spike_times = [spikes_x,spikes_t];
        %Find Events
        events = []; 
        event_lengths = [];
        last_start = spikes_t(1);
        last_time = spikes_t(1);
        spike_count = 1;
        for t_i = 2:length(spikes_t)
            s_i = spikes_t(t_i);
            if s_i - last_time <= parameters.IES
                last_time = s_i;
                spike_count = spike_count + 1;
            else
                if (last_start ~= last_time) && (spike_count > parameters.event_cutoff*parameters.n) %weed out events w/ too few spikes
                    events = [events; [last_start, last_time]]; %#ok<AGROW> %add the last range of spikes to the events vector
                    event_lengths = [event_lengths, (last_time - last_start)*parameters.dt]; %#ok<*AGROW>
                end
                last_start = s_i;
                last_time = s_i;
                spike_count = 1;
            end
        end
        if (last_start ~= last_time) && (spike_count > parameters.event_cutoff*parameters.n) %weed out events w/ too few spikes
            events = [events; [last_start, last_time]]; %#ok<AGROW> %add the last interval
            event_lengths = [event_lengths, (last_time - last_start)*parameters.dt];
        end
        clear t_i s_i last_start last_time
        [num_events,~] = size(events);
        %Store Event Info
        network_spike_sequences_no_I(i).num_events = num_events;
        network_spike_sequences_no_I(i).event_lengths = event_lengths;
        network_spike_sequences_no_I(i).events = events;
        network_cluster_sequences_no_I(i).num_events = num_events;
        network_cluster_sequences_no_I(i).event_lengths = event_lengths;
        network_cluster_sequences_no_I(i).events = events;
        %If events are not empty, store spike order, ranks, and nonspiking
        %neurons
        if ~isempty(events)
            for j = 1:num_events %Run through each event
                start_ind = find(spike_times(:,2) == events(j,1));
                end_ind = find(spike_times(:,2) == events(j,2));
                event_spike_neurons = spike_times(start_ind:end_ind,1);
                neuron_order = unique(event_spike_neurons,'stable');
                ranks = zeros(length(e_indices),1);
                for k = 1:length(e_indices)
                    rank_val = find(neuron_order == k,1);
                    if ~isempty(rank_val)
                        ranks(k) = rank_val;
                    end
                end    
                clear k rank_val
                nonspiking_neurons = isnan(ranks./ranks);
                network_spike_sequences_no_I(i).spike_order.(strcat('sequences_',string(j))) = neuron_order;
                network_spike_sequences_no_I(i).spike_ranks.(strcat('sequences_',string(j))) = ranks;
                network_spike_sequences_no_I(i).nonspiking_neurons.(strcat('sequences_',string(j))) = nonspiking_neurons;
            end
        end
        clear j start_ind end_ind event_spike_neurons neuron_order ranks ...
            k rank_val nonspiking_neurons
        %Separately store cluster sequence information
        bin_width = 5*10^(-3); %5 ms bin
        bin_size = ceil(bin_width/parameters.dt); %number of timesteps to use in a bin
        for e_i = 1:num_events
            cluster_spikes = network.cluster_mat(:,e_indices)*spikes_V_m(:,events(e_i,1):events(e_i,2));
            cluster_mov_sum = movsum(cluster_spikes',bin_size)';
            normalized_cluster_spikes = cluster_spikes ./ sum(cluster_spikes,1);
            normalized_cluster_spikes(isnan(normalized_cluster_spikes)) = 0;
            normalized_cluster_mov_sum = cluster_mov_sum ./ sum(cluster_mov_sum,1);
            normalized_cluster_mov_sum(isnan(normalized_cluster_mov_sum)) = 0;
            network_cluster_sequences_no_I(i).clusters.(strcat('sequence_',string(e_i))) = cluster_spikes;
            network_cluster_sequences_no_I(i).movsum.(strcat('sequence_',string(e_i))) = cluster_mov_sum;
            network_cluster_sequences_no_I(i).normalized_clusters.(strcat('sequence_',string(e_i))) = normalized_cluster_spikes;
            network_cluster_sequences_no_I(i).normalized_cluster_mov_sum.(strcat('sequence_',string(e_i))) = normalized_cluster_mov_sum;
        end
        clear bin_width bin_size cluster_spikes cluster_mov_sum ...
            normalized_cluster_spikes normalized_cluster_mov_sum e_i
    end
end