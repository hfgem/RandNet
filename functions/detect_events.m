function [network_spike_sequences, network_cluster_sequences, overallResults] = detect_events(parameters, ...
                    network, V_m , ithTest, network_spike_sequences, network_cluster_sequences, varargin)
    %_________
    %ABOUT: This function detects events from the activity in V_m and the
    %event criterion parameters in "paramters"
    %
    %INPUTS:
    %   parameters = simulations' parameter structure
    %   network = newtork structure
    %   V_m = a single simulations membrane potential matrix
    %   ithTest = if network is being simulated multiple times, this saves
    %   results to the ith component of the output structures
    %   network_spike_sequences, network_cluster_sequences = initialized or
    %   pre-existing structures that have results for the ithTest appended
    %
    %OUTPUTS:
    %   network_spike_sequences = struct of spike sequences
    %   network_cluster_sequences = struct of cluster sequences
    %_________
    
    % Defaults for optional parameters
    E_events_only = 0;   % if 1, only consider E_cells for events
    
    % Read in optional parameters, to overwrite above defaults
    for i=1:2:length(varargin)
        switch varargin{i}
            case 'E_events_only'
                E_events_only = varargin{i+1};
            otherwise
                error('detect_events: Unknown input')
        end
    end
    
    %Find spike profile
    if E_events_only
        indices = network.E_indices;
    else
        indices = network.all_indices;
    end
    spikes_V_m = V_m(indices,:) >= parameters.V_th;
    
    [spikes_x,spikes_t] = find(spikes_V_m);
    % max_time = max(spikes_t);
    spiking_neurons = unique(spikes_x, 'stable');
    
    avg_event_length = nan; % Initialize, in case there are no events
    event_lengths = [];
    events = [];

    %Find maximum firing rate + average maximum firing rates of neurons
    all_fr = sum(spikes_V_m,2)/parameters.t_max;
    max_fr = max(all_fr);
    avg_fr = mean(all_fr);
    %display(avg_fr)
        
    %TEST 1: The number of neurons participating in a sequence must
    %pass a threshold:
    if length(spiking_neurons) >= parameters.event_cutoff*numel(indices)

        %TEST 2: The firing rate must fall within a realistic range
        if and(avg_fr>= parameters.min_avg_fr, avg_fr <= parameters.max_avg_fr)
            %Find event times
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
                    if (last_start ~= last_time) && (spike_count > parameters.event_cutoff*numel(indices)) %weed out events w/ too few spikes
                        events = [events; [last_start, last_time]]; %#ok<AGROW> %add the last range of spikes to the events vector
                        event_lengths = [event_lengths, (last_time - last_start)*parameters.dt]; %#ok<*AGROW>
                    end
                    last_start = s_i;
                    last_time = s_i;
                    spike_count = 1;
                end
            end
            if (last_start ~= last_time) && (spike_count > parameters.event_cutoff*numel(indices)) %weed out events w/ too few spikes
                events = [events; [last_start, last_time]]; %#ok<AGROW> %add the last interval
                event_lengths = [event_lengths, (last_time - last_start)*parameters.dt]; %#ok<*SAGROW>
            end
            avg_event_length = mean(event_lengths);
            [num_events,~] = size(events);
            %save to both structures
            network_spike_sequences(ithTest).events = events;
            network_spike_sequences(ithTest).event_lengths = event_lengths;
            network_cluster_sequences(ithTest).events = events;
            network_cluster_sequences(ithTest).event_lengths = event_lengths;

            %TEST 3: The sequence(s) of firing is(are) within
            %reasonable lengths.
            if and(avg_event_length >= parameters.min_avg_length, avg_event_length <= parameters.max_avg_length)

                %Find spike sequences
                for e_i = 1:num_events
                    %store spike orders for each event
                    event_spikes = spikes_V_m(:,events(e_i,1):events(e_i,2));
                    [e_spikes_x, ~] = find(event_spikes);
                    spike_order = unique(e_spikes_x,'stable');
                    network_spike_sequences(ithTest).spike_order.(strcat('sequence_',string(e_i))) = spike_order;
                    %store ranks for each neuron
                    ranks_vec = nan(1,numel(indices));
                    for k = 1:length(spike_order)
                        n_ind = spike_order(k);
                        ranks_vec(1,n_ind) = k;
                    end
                    network_spike_sequences(ithTest).spike_ranks.(strcat('sequence_',string(e_i))) = ranks_vec;
                    %store nonspiking neurons
                    nonspiking_neurons = isnan(ranks_vec./ranks_vec);
                    network_spike_sequences(ithTest).nonspiking_neurons.(strcat('sequence_',string(e_i))) = nonspiking_neurons;
                end
                clear e_i event_spikes e_spikes_x spike_order ranks_vec k n_ind nonspiking_neurons


                %Find cluster sequence per event by moving bin
                bin_size = ceil(parameters.bin_width/parameters.dt); %number of timesteps to use in a bin
                for e_i = 1:num_events
                    cluster_spikes = network.cluster_mat(:,indices)*spikes_V_m(:,events(e_i,1):events(e_i,2));
                    cluster_mov_sum = movsum(cluster_spikes',bin_size)';
                    network_cluster_sequences(ithTest).clusters.(strcat('sequence_',string(e_i))) = cluster_spikes;
                    network_cluster_sequences(ithTest).movsum.(strcat('sequence_',string(e_i))) = cluster_mov_sum;
                end
                clear bin_width bin_size cluster_spikes cluster_mov_sum e_i
            end %Sequence length loop
        end %Firing rate loop
    end %Number of neurons in sequence loop
    
    % Statistics over the whole simulation
    overallResults(1) = length(spiking_neurons); % number of spiking neurons
    overallResults(2) = avg_fr; %average firing rate
    overallResults(3) = avg_event_length; %Average event length
    overallResults(4) = size(events, 1); % Number of events
    
    if size(events, 1)==0 % set results as empty, to prevent errors upon analysis
        network_spike_sequences(ithTest).events = [];
        network_spike_sequences(ithTest).event_lengths = [];
        network_spike_sequences(ithTest).spike_order = [];
        network_spike_sequences(ithTest).ranks_vec = [];
        network_spike_sequences(ithTest).nonspiking_neurons = [];
    end
end