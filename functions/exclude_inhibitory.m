%Generate sequence file excluding inhibitory neurons.

function sequences = exclude_inhibitory(original_sequences, e_indices, parameters)
    %Pull Spikes
    sequences = struct;
    for i = 1:length(original_sequences)
        %Find and store all spike times, excluding inhibitory neurons
        V_m = original_sequences(i).V_m;
        spikes_V_m = V_m(e_indices,:) >= parameters.V_th;
        sequences(i).spikes_V_m = spikes_V_m;
        %Find events
        [spikes_x,spikes_t] = find(spikes_V_m);
        spike_times = [spikes_x,spikes_t];
        spiking_neurons = unique(spikes_x, 'stable');
        sequences(i).spiking_neurons = spiking_neurons;
        sequences(i).spike_times = spike_times;
        %Find Events
        events = []; 
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
        end
        [num_events,~] = size(events);
        sequences(i).num_events = num_events;
        sequences(i).events = events;
        %If events are not empty, continue
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
                sequences(i).ranks.(strcat('sequences_',string(j))) = ranks;
                sequences(i).neuron_order.(strcat('sequences_',string(j))) = neuron_order;
            end
        end
    end
end