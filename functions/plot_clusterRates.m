function plot_clusterRates(spikes_V_m, parameters, network, varargin)
% Generates a plot of the mean spike rate, mean cluster-wise spike rate, 
% and z-scored cluster-wise spike rate from the data in spikes_V_m
% 
% spikes_V_m: binary spike matrix from one trial
% parameters: the parameter structure for the simulation
% network: the network structure
%
% Example usage, after running a simulation from randnet.m:
% smoothWindow = 10 * (1/parameters.dt * 1/1000); %gaussian kernel width for smoothing firing rate curves
% plot_clusterRates(spikes_V_m, parameters, network, 'smoothWindow', smoothWindow);


%% Default parameters:
E_only = 1; % if 1, only include E-cells
legend_flag = 0; % if 1, include cluster legend
smoothWindow = 50 * (1/parameters.dt * 1/1000); %gaussian kernel width for smoothing firing rate curves

%% Read in optional parameters, to overwrite above defaults
for i=1:2:length(varargin)
    switch varargin{i}
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

t = [0:parameters.dt:parameters.t_max];


if E_only
    spike_mat = spikes_V_m(network.E_indices,:);
else
    spike_mat = spikes_V_m;
end

y = zeros(parameters.clusters, size(spike_mat, 2)); % num spikes each cluster fired each time step
for ithCluster = 1:parameters.clusters
    clusterMember = network.cluster_mat(ithCluster,network.E_indices);
    y(ithCluster,:) = clusterMember*spike_mat/sum(clusterMember)/parameters.dt;
end

figure; 

yRate = smoothdata(mean(spike_mat, 1)/parameters.dt, 'gaussian', smoothWindow);
ax1 = subplot(3, 1, 1); plot(t, yRate);
ylabel('Pop. mean rate (Hz)')

ySmoothed = smoothdata(y, 2, 'gaussian', smoothWindow);
ax2 = subplot(3, 1, 2); plot(t, ySmoothed)
ylabel('Cluster mean rate (Hz)')
if legend_flag
    legend( sprintfc('Cluster %g', 1:parameters.clusters), 'Location', 'Best' )
end

yZScore = (ySmoothed-mean(ySmoothed, 2))./std(ySmoothed, [], 2) ; % z-score across time, for each cluster
% yZScore = (ySmoothed-mean(ySmoothed, 1))./std(ySmoothed, [], 1) ; % z-score across clusters, for each time point
ax3 = subplot(3, 1, 3); plot(t, yZScore);;
ylabel('Cluster z-score')
xlabel('Time (s)'); 
if legend_flag
    legend( sprintfc('Cluster %g', 1:parameters.clusters), 'Location', 'Best' )
end
    

linkaxes([ax1 ax2 ax3],'x')

end