function plot_sequence_clusterRates(network_spike_sequences, plotIDs, spikes_V_m, parameters, network, varargin)
% For each event in plotIDs, it plots the event raster, the cluster-wise
% mean spike rate, and the cluster-wise z-scored rate
% 
% network_spike_sequences: the structure created by detect_PBE
% plotIDs: vector of the sequence IDs to plot
% spikes_V_m: binary spike matrix from one trial
% parameters: the parameter structure for the simulation
% network: the network structure
%
% Example usage, after running a simulation from randnet.m:
% [network_spike_sequences] = detect_PBE(spikes_V_m, parameters);
% plotIDs = [1, 3, 4]; % indexes of the sequences to plot
% plot_sequence_clusterRates(network_spike_sequences, plotIDs, spikes_V_m, parameters, network);


%% Default parameters:
ithTrial = 1;
E_only = 1; % if 1, only include E-cells
legend_flag = 0; % if 1, include cluster legend
smoothWindow = 5 * (1/parameters.dt * 1/1000); %gaussian kernel width for smoothing firing rate curves

%% Read in optional parameters, to overwrite above defaults
for i=1:2:length(varargin)
    switch varargin{i}
        case 'ithTrial'
            ithTrial = varargin{i+1};
        case 'E_only'
            E_only = varargin{i+1};
        case 'legend_flag'
            legend_flag = varargin{i+1};
        case 'smoothWindow'
            smoothWindow = varargin{i+1};
        otherwise
            error('plot_sequence_clusterRates: Unknown input')
    end
end

%%

events = network_spike_sequences.events;
num_events = size(events, 1);

assert([num_events>=max(plotIDs)], 'largest plotIDs is larger than the number of events')

for e_i = plotIDs
    
    spike_ranks = network_spike_sequences(ithTrial).ranks_vec(:,e_i);
    [~, spike_order] = sort(spike_ranks);
    [~, Ie] = sort(spike_ranks(network.E_indices));
    [~, Ii] = sort(spike_ranks(network.I_indices));
    
    event_spikes = spikes_V_m(:,events(e_i,1):events(e_i,2));
    
    if ~E_only
        reordered_spikes = event_spikes(spike_order,:);
    else
        reordered_spikes = [spikes_V_m(network.E_indices(Ie),events(e_i,1):events(e_i,2)); ...
                        spikes_V_m(network.I_indices(Ii),events(e_i,1):events(e_i,2))];
    end
            
    figure;
    sgtitle(['Event ' num2str(e_i)])
    subplot(1, 3, 1)
    plotSpikeRaster( reordered_spikes, 'TimePerBin', parameters.dt*1000, 'PlotType', 'scatter');
    spikeAx = gca;
    xlabel('Time (ms)'); ylabel('Neurons, sorted')

    
    y = zeros(parameters.clusters, size(event_spikes, 2)); % num spikes each cluster fired each time step
    for iCluster = 1:parameters.clusters
        clusterMember = network.cluster_mat(iCluster,network.E_indices);
        eCellSpikes = event_spikes(network.E_indices,:);
        y(iCluster,:) = clusterMember*eCellSpikes;
    end
    
    ySmoothed = smoothdata(y, 2, 'gaussian', smoothWindow);
  
    subplot(1, 3, 2)
    plot([1:size(event_spikes, 2)]*parameters.dt*1000, ySmoothed)
    xlabel('Time (ms)'); ylabel('Cluster firing rate')
    xlim(spikeAx.XLim)
    if legend_flag
        legend( sprintfc('Cluster %g', 1:parameters.clusters), 'Location', 'Best' )
    end
    
    yZScore = (ySmoothed-mean(ySmoothed, 2))./std(ySmoothed, [], 2) ;
    subplot(1, 3, 3)
    plot([1:size(event_spikes, 2)]*parameters.dt*1000, yZScore)
    xlabel('Time (ms)'); ylabel('Cluster Z-score')
    xlim(spikeAx.XLim)
    if legend_flag
        legend( sprintfc('Cluster %g', 1:parameters.clusters), 'Location', 'Best' )
    end
    
end

end