%% Plot various analyses from results of RandNet simulation
%

seed = 1;
ithTrial = 1;

E_spikes_V_m = V_m(network.E_indices,:) >= parameters.V_th;
spikes_V_m = V_m >= parameters.V_th;
t = [0:parameters.dt:parameters.t_max];

rng(1)
PFpeaksSequence = randperm(parameters.n_E)'; % dummy data, for testing

%% Main plots
plot_randnet_results(parameters, network, V_m, G_in, network_spike_sequences, ithTest, net_save_path)


%% Plot multiple sequences with different sortings
[network_spike_sequences] = detect_PBE(E_spikes_V_m, parameters);
plotIDs = [1, 2]; % indexes of the sequences to plot
plot_ixj_sequences(network_spike_sequences, plotIDs, spikes_V_m, parameters, network);


%% Use dimensionality reduction techniques to plot network structure
E_only = 1; % only plot E-E connections
dimRedInput = 'W'; % 'W', 'WW', 'normClust', 'clust' 
scatterSize = [1, 1];
plot_dimRedNet(network, dimRedInput, 'seed', seed, 'E_only', E_only, 'scatterSize', scatterSize);


%% Plot cluster-wise activation, for each event
[network_spike_sequences] = detect_PBE(E_spikes_V_m, parameters);
plotIDs = [1, 2]; % indexes of the sequences to plot
plot_sequence_clusterRates(network_spike_sequences, plotIDs, spikes_V_m, parameters, network);


%% Plot cluster-wise spike rate over entire simulation
smoothWindow = 100 * (1/parameters.dt * 1/1000); %gaussian kernel width for smoothing firing rate curves
plot_clusterRates(E_spikes_V_m, parameters, network, 'smoothWindow', smoothWindow);


%% Plots sequence-sequence correlation analysis
seed = 1;
ithTrial = 1;
correlationType = 'Pearson'; % Pearson, Kendall, Spearman
usRelRank = 1;
maxNClust = 4;

[network_spike_sequences] = detect_PBE(E_spikes_V_m, parameters);
ranks_vec = network_spike_sequences.ranks_vec;
plot_seq_seq_corr(ranks_vec, 'correlationType', correlationType);


%% Matching index for sequence-sequence comparisons
seed = randi(10^6);
nShuffles= 100;
maxNClust = 5;
penalize_nonspike = 0;

[network_spike_sequences] = detect_PBE(E_spikes_V_m, parameters);
ranks_vec = network_spike_sequences(ithTrial).ranks_vec;
plot_seq_seq_MI(ranks_vec)


%% Pearson correlation for sequence-PF comparisons
% Requires place fields
if exist('PFpeaksSequence')
    correlationType = 'Pearson'; % Pearson, Kendall, Spearman
    [network_spike_sequences] = detect_PBE(E_spikes_V_m, parameters);
    ranks_vec = network_spike_sequences(ithTrial).ranks_vec; % ranks for each detected PBE
    plot_PF_seq_corr(ranks_vec, PFpeaksSequence, 'correlationType', correlationType);
end


%% Mean relative rank vs place field rank
% Requires place fields
if exist('PFpeaksSequence')
    minPartic = 1; % minimum number of sequences a cell needs to participate in
    [network_spike_sequences] = detect_PBE(E_spikes_V_m, parameters);
    ranks_vec = network_spike_sequences(ithTrial).ranks_vec;
    plot_meanRelRank(ranks_vec, PFpeaksSequence)
end


%% Pairwise spike probabilities
% Requires place fields
if exist('PFpeaksSequence')
    plot_pairwise_spikeProb(E_spikes_V_m, PFpeaksSequence)
end
