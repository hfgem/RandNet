function [network_spike_sequences] = detect_PBE(spikes, parameters, varargin)
    %_________
    %ABOUT: Detects population burst events based on criterion in the
    %parameters structure
    %
    %INPUTS:
    %   spikes = nNeurons x nTimesteps binary spike matrix
    %   parameters = simulation and analysis parameter structure
    %
    %OUTPUTS:
    %   network_spike_sequences = struct of spike sequences
    %       Fields: 
    %_________

    
    %% Read in optional parameters, to overwrite above defaults
    for i=1:2:length(varargin)
        switch varargin{i}
            otherwise
                error('detect_events: Unknown input')
        end
    end
       
    %% Detect PBE timepoints
    t = 0:parameters.dt:(size(spikes, 2)-1)*parameters.dt;
    meanPopRate = smoothdata( mean(spikes, 1)/parameters.dt, 'gaussian', parameters.PBE_window);

    PBEthresh = max(mean(meanPopRate)+(parameters.PBE_zscore*std(meanPopRate)), parameters.PBE_min_Hz); % threshold rate for PBE detection
    % figure; plot(t, meanPopRate); hold on; yline(PBEthresh)
    % figure; plotSpikeRaster(spikes, 'TimePerBin', parameters.dt, 'PlotType', 'scatter'); 
    
    
    PBE_candidate = [0, meanPopRate>PBEthresh]; % binary vector of candidate PBE times
    % figure; plot(t, PBE_candidate)
    
    % combine candidate PBE that are close in time
    onsetInds = find(diff(PBE_candidate)==1)+1; % indexes where a candidate PBE starts
    offsetInds = find(diff(PBE_candidate)==-1)+1; % indexes where a candidate PBE ends
    for i = numel(offsetInds)-1:-1:1
        if onsetInds(i+1)-offsetInds(i)<parameters.PBE_max_combine/parameters.dt
            PBE_candidate(offsetInds(i)-1:onsetInds(i+1)+1) = 1;
        end
    end
    % figure; plot(t, PBE_candidate)
    
    % Remove candidate PBEs that don't last long enough
    onsetInds = find(diff(PBE_candidate)==1)+1; % indexes where a candidate PBE starts
    for i = numel(onsetInds):-1:1
        if onsetInds(i)+parameters.PBE_min_dur/parameters.dt > numel(t)
            PBE_candidate(onsetInds(i):end) = 0;
        elseif ~all(PBE_candidate(onsetInds(i):onsetInds(i)+parameters.PBE_min_dur/parameters.dt))
            if i~=numel(onsetInds)
                PBE_candidate(onsetInds(i):onsetInds(i+1)-1) = 0;
            else
                PBE_candidate(onsetInds(i):onsetInds(i)+parameters.PBE_min_dur/parameters.dt) = 0;
            end
        end
    end
    % figure; plot(t, PBE_candidate)    
    
    
    
    %% Analyze PBE contents
    
    onsetInds_final = find(diff(PBE_candidate)==1)+1; % indexes where processed PBEs starts
    offsetInds_final = find(diff(PBE_candidate)==-1)+1; % indexes where processed PBEs ends
    
    if numel(onsetInds_final) > numel(offsetInds_final) % if event continues to the end of the trial, terminate it then
        offsetInds_final = [offsetInds_final, numel(meanPopRate)];
    end
    
    events = [onsetInds_final', offsetInds_final']; 
    event_lengths = [offsetInds_final' - onsetInds_final']*parameters.dt; 
    
    if size(events, 1)==0 % If no events, set results as empty
        network_spike_sequences.events = [];
        network_spike_sequences.event_lengths = [];
        network_spike_sequences.spike_order = {};
        network_spike_sequences.ranks_vec = [];
        network_spike_sequences.frac_spike = {};
        
    else % If PBEs detected, calculate and store results
        
        network_spike_sequences.events = events;
        network_spike_sequences.event_lengths = event_lengths;
        
        %Find spike sequence for each event
        for ithEvent = 1:size(events, 1)

            event_spikes = spikes(:,events(ithEvent,1):events(ithEvent,2));
            
            if parameters.plotResults == 1 %Plot event raster
                figure;
                plotSpikeRaster(event_spikes, 'TimePerBin', parameters.dt, 'PlotType', 'scatter'); 
            end
            
            [e_spikes_x, ~] = find(event_spikes); % find neuron index of all spikes in order of spike-time
            spike_order = unique(e_spikes_x,'stable'); % take only first spike of each neuron
            network_spike_sequences.spike_order{ithEvent} = spike_order;
            
            %store ranks for each neuron
            ranks_vec = nan(1, size(spikes, 1));
            for k = 1:length(spike_order)
                n_ind = spike_order(k);
                ranks_vec(1,n_ind) = k;
            end
            network_spike_sequences.ranks_vec(:,ithEvent) = ranks_vec;
            network_spike_sequences.frac_spike{ithEvent} = sum(~isnan(ranks_vec))/numel(ranks_vec);
        end
                
    end
        
end