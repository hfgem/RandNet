function [network_spike_sequences, overallResults] = detect_events(parameters, ...
                    spikes_V_m, ithParamSet, ithNet, ithTest)
    %_________
    %ABOUT: This function detects events from the activity in V_m and the
    %event criterion parameters in "parameters"
    %
    %INPUTS:
    %   parameters = simulations' parameter structure
    %   spikes_V_m = a single simulation's binary spikes matrix - can be
    %       passed in limited to any subpopulation of interest
    %   ind_test = if the results are to be plotted or saved, this will set
    %       the save name
    %OUTPUTS:
    %   network_spike_sequences = struct of spike sequences
    %   network_cluster_sequences = struct of cluster sequences
    %_________
    
    %Find spike profile
    [num_neur,~] = size(spikes_V_m);
    [spikes_x,spikes_t] = find(spikes_V_m);
    spiking_neurons = unique(spikes_x, 'stable');
    
    event_lengths = [];
    events = [];

    %Find maximum firing rate + average maximum firing rates of neurons
    all_fr = sum(spikes_V_m,2)/parameters.t_max;
    avg_fr = mean(all_fr);
    %display(avg_fr)
    
    network_spike_sequences = struct;
    
    if length(spikes_t) > 1
        %Find event times
        last_start = spikes_t(1);
        last_time = spikes_t(1);
        spike_count = 1;
        for t_i = 2:length(spikes_t)
            s_i = spikes_t(t_i);
            if s_i - last_time <= parameters.IES %check if the next spike is within an event interval
                last_time = s_i; %store the updated last time of a spike in the event
                spike_count = spike_count + 1; %update the count of spikes in the event
            else %if the spike is not within an event timebin update parameters to store a new event
                if (last_start ~= last_time) && (spike_count >= parameters.event_cutoff*num_neur) %weed out events w/ too few spikes
                    events = [events; [last_start, last_time]]; %#ok<AGROW> %add the last range of spikes to the events vector
                    event_lengths = [event_lengths, (last_time - last_start)*parameters.dt]; %#ok<*AGROW>
                end
                last_start = s_i;
                last_time = s_i;
                spike_count = 1;
            end
        end
        %Add last event if it has enough neurons
        if (last_start ~= last_time) && (spike_count >= parameters.event_cutoff*num_neur)
            events = [events; [last_start, last_time]]; %#ok<AGROW> %add the last interval
            event_lengths = [event_lengths, (last_time - last_start)*parameters.dt]; %#ok<*SAGROW>
        end

        %Calculate event statistics
        avg_event_length = mean(event_lengths);
        [num_events,~] = size(events);

        % If no events, set results as empty
        if size(events, 1)==0
            network_spike_sequences.events = [];
            network_spike_sequences.event_lengths = [];
            network_spike_sequences.spike_order = {};
            network_spike_sequences.frac_spike = {};
            network_spike_sequences.ranks_vec = [];
        else
            %save to both structures
            network_spike_sequences.events = events;
            network_spike_sequences.event_lengths = event_lengths;

            if parameters.plotResults == 1
                f = figure(); %#ok<NASGU>
            end
            
            any_success = 0;

            %Find spike sequence for each event
            for ithEvent = 1:num_events
                %store spike orders for each event
                event_spikes = spikes_V_m(:,events(ithEvent,1):events(ithEvent,2));

                if parameters.plotResults == 1 %Plot event raster
                    subplot(1,num_events,ithEvent)
                    plotSpikeRaster(event_spikes, 'TimePerBin', parameters.dt, 'PlotType', 'scatter');
                    title(sprintf('Event %i Raster',ithEvent))
                end

                [e_spikes_x, ~] = find(event_spikes);
                spike_order = unique(e_spikes_x,'stable');
                network_spike_sequences.spike_order{ithEvent} = spike_order;

                %store ranks for each neuron
                ranks_vec = nan(1,num_neur);
                for k = 1:length(spike_order)
                    n_ind = spike_order(k);
                    ranks_vec(1,n_ind) = k;
                end
                network_spike_sequences.ranks_vec(:,ithEvent) = ranks_vec;
                network_spike_sequences.frac_spike{ithEvent} = sum(~isnan(ranks_vec))/num_neur;

                partic_check = length(spike_order) >= parameters.event_cutoff;
                fr_check = parameters.min_avg_fr <= avg_fr <= parameters.max_avg_fr;
                len_check = parameters.min_avg_length <= avg_event_length <= parameters.max_avg_length;
                overall_check = partic_check*fr_check*len_check;
                if overall_check ~= 0
                    any_success = 1;
                end
                
                %plot only  if meets criteria
                if parameters.plotResults == 1
                    if overall_check ~= 0
                        subplot(1,num_events,ithEvent)
                        plotSpikeRaster(event_spikes, 'TimePerBin', parameters.dt, 'PlotType', 'scatter');
                        title(sprintf('Event %i Raster',ithEvent))
                    else
                        subplot(1,num_events,ithEvent)
                        title(sprintf('Event %i Does Not Meet Criteria',ithEvent))
                    end
                end
            end    

            if (any_success == 1) && (parameters.plotResults == 1) && (parameters.saveFlag == 1) %Save the raster plots of events
                f.Units = 'Normalized';
                f.OuterPosition = [0 0 1 1];
                if ~isfolder(strcat(parameters.save_path,'/sim_plots'))
                    mkdir(parameters.save_path,'/sim_plots')
                end
                savefig(f,strcat(parameters.save_path,'/sim_plots','/sim_',string(ithParamSet),'_',string(ithNet),'_',string(ithTest),'_events.fig'))
                saveas(f,strcat(parameters.save_path,'/sim_plots','/sim_',string(ithParamSet),'_',string(ithNet),'_',string(ithTest),'_events.jpg'))
                close(f)
            elseif (any_success == 0) && (parameters.plotResults == 1)
                close(f)
            end
        end
        
        % Statistics over the whole simulation
        overallResults(1) = length(spiking_neurons); % number of spiking neurons
        overallResults(2) = avg_fr; %average firing rate
        overallResults(3) = avg_event_length; %Average event length
        overallResults(4) = size(events, 1); % Number of events

    else
        network_spike_sequences.events = [];
        network_spike_sequences.event_lengths = [];
        network_spike_sequences.spike_order = {};
        network_spike_sequences.frac_spike = {};
        network_spike_sequences.ranks_vec = [];
        
        % Statistics over the whole simulation
        overallResults(1) = length(spiking_neurons); % number of spiking neurons
        overallResults(2) = avg_fr; %average firing rate
        overallResults(3) = 0; %Average event length
        overallResults(4) = size(events, 1); % Number of events

    end %at least 1 spike if statement
    
end