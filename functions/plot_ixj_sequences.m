function figHandle = plot_ixj_sequences(network_spike_sequences, plotIDs, spikes_V_m, parameters, network, varargin)
% Plots plotIDs sequences from network_spike_sequences(ithTrial).spike_ranks
% sorted against each of the other plotIDs sequences in an nToPlot by nToPlot
% grid of subplots
% 
% network_spike_sequences: the structure created by detect_PBE
% plotIDs: vector of the sequence IDs to plot
% spikes_V_m: binary spike matrix from one trial
% parameters: the parameter structure for the simulation
% ithTrial: specifices the index of network_spike_sequences(ithTrial) to use
%     if spikes_V_m does not correspond to network_spike_sequences(1).
%     Optional input parameter.
%
% figHandle: optional output, figure handle of the generated figure
%
% Example usage, after running a simulation from randnet.m:
% [network_spike_sequences] = detect_PBE(spikes_V_m, parameters); % detect PBEs
% plotIDs = [1, 2, 4]; % indexes of the sequences to plot
% plot_ixj_sequences(network_spike_sequences, plotIDs, spikes_V_m, parameters); % generate figure


%% Default parameters:
ithTrial = 1;
separate_EI = 1; % if 1 , separate E and I cells on raster

%% Read in optional parameters, to overwrite above defaults
for i=1:2:length(varargin)
    switch varargin{i}
        case 'ithTrial'
            ithTrial = varargin{i+1};
        case 'separate_EI'
            separate_EI = varargin{i+1};
        otherwise
            error('plot_ixj_sequences: Unknown input')
    end
end

%%

events = network_spike_sequences.events;
num_events = size(events, 1);

assert([num_events>=max(plotIDs)], 'largest plotIDs is larger than the number of events')

figHandle = figure;
ithPlot = 1; % initialize to 1, then increment inside loop

for k = plotIDs
    
    spike_ranks = network_spike_sequences(ithTrial).ranks_vec(:,k);
    [~, spike_order] = sort(spike_ranks);
    [~, Ie] = sort(spike_ranks(network.E_indices));
    [~, Ii] = sort(spike_ranks(network.I_indices));

    for e_i = plotIDs
        
        event_spikes = spikes_V_m(:,events(e_i,1):events(e_i,2));
        
        if ~separate_EI
            reordered_spikes = event_spikes(spike_order,:);
        else
            reordered_spikes = [spikes_V_m(network.E_indices(Ie),events(e_i,1):events(e_i,2)); ...
                            spikes_V_m(network.I_indices(Ii),events(e_i,1):events(e_i,2))];
        end
        
        subplot(numel(plotIDs),numel(plotIDs),ithPlot); 
        plotSpikeRaster( reordered_spikes, 'TimePerBin', parameters.dt*1000, 'PlotType', 'scatter');
        
        if ismember(ithPlot, [7, 8, 9])
            xlabel('Time (ms)')
        end
        if e_i==1
            ylabel( [{'Cells sorted'}, {strcat('by Event #',string(k))}] )
        end
        if k==1
            title( strcat('Event #',string(e_i)) )
        end
        
        ithPlot = ithPlot + 1;
    end
    sgtitle(['Plotting events ', strjoin(cellstr(num2str(plotIDs')),','), ' of ', ...
        num2str(num_events), ' total events'])
    
end

end