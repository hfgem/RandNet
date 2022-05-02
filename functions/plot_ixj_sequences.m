function figHandle = plot_ixj_sequences(network_spike_sequences, plotIDs, spikes_V_m, parameters, network, varargin)
% Plots plotIDs' sequences from network_spike_sequences(ithTrial).spike_ranks
% sorted against each of the other plotIDs sequences in an nToPlot by nToPlot
% grid of subplots
% 
% network_spike_sequences: the structure created by detect_PBE
% plotIDs: vector of the sequence IDs to plot
% spikes_V_m: binary spike matrix from one trial. Can be either E-cell only
%     or E and I cells together, even if detect_PBE network_spike_sequences 
%     is only E-cells.
% parameters: the parameter structure for the simulation
% ithTrial: specifices the index of network_spike_sequences(ithTrial) to use
%     if spikes_V_m does not correspond to network_spike_sequences(1).
%     Optional input parameter.
%
% figHandle: optional output, figure handle of the generated figure
%
% Example usage, after running a simulation from randnet.m:
% [network_spike_sequences] = detect_PBE(E_spikes_V_m, parameters); % detect PBEs based on E-cells only
% plotIDs = [1, 2, 4]; % indices of the sequences to plot
% plot_ixj_sequences(network_spike_sequences, plotIDs, spikes_V_m, parameters); % generate figure, with E-cells sorted and I-cells unsorted


%% Default parameters:
ithTrial = 1;

%% Read in optional parameters, to overwrite above defaults
for i=1:2:length(varargin)
    switch varargin{i}
        case 'ithTrial'
            ithTrial = varargin{i+1};
        otherwise
            error('plot_ixj_sequences: Unknown input')
    end
end

%% Main function:

E_only = [size(network_spike_sequences.ranks_vec,1)==parameters.n_E]; % if true, detect_PBE only used E-cells

events = network_spike_sequences.events;
num_events = size(events, 1);

assert([num_events>=max(plotIDs)], 'largest plotIDs is larger than the number of events')

figHandle = figure;
ithSubPlot = 1; % initialize to 1, then increment inside loop
for k = plotIDs
    
    spike_ranks = network_spike_sequences(ithTrial).ranks_vec(:,k);
    
    if E_only 
        [~, Ie] = sort(spike_ranks);

        if [size(spikes_V_m,1)==size(network_spike_sequences.ranks_vec,1)]
            I_outer = Ie; 
        else
            I_outer = [network.E_indices(Ie), network.I_indices];
        end
    else
        [~, Ie] = sort(spike_ranks(network.E_indices));
        [~, Ii] = sort(spike_ranks(network.I_indices));
        I_outer = [network.E_indices(Ie), network.I_indices(Ii)];
    end

    for e_i = plotIDs
        
        reordered_spikes = [spikes_V_m(I_outer,events(e_i,1):events(e_i,2))];
        
        subplot(numel(plotIDs),numel(plotIDs),ithSubPlot); 
        plotSpikeRaster( reordered_spikes, 'TimePerBin', parameters.dt*1000, 'PlotType', 'scatter');
        
        % Plot labels, conditional on subplot
        if ismember(ithSubPlot, [7, 8, 9])
            xlabel('Time (ms)')
        end
        if e_i==1
            ylabel( [{'Cells sorted'}, {strcat('by Event #',string(k))}] )
        end
        if k==1
            title( strcat('Event #',string(e_i)) )
        end
        
        ithSubPlot = ithSubPlot + 1;
    end
    sgtitle(['Plotting events ', strjoin(cellstr(num2str(plotIDs')),','), ' of ', ...
        num2str(num_events), ' total events'])
    
end

end