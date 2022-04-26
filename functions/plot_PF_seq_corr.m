function plot_PF_seq_corr(ranks_vec, PFpeaksSequence, varargin)
% Calculates correlation of sequences in network_spike_sequences to the
% Place field sequence in PFpeaksSequence
% Also calculates shuffle sequence correlations to PFpeakSequence
% 
% ranks_vec: the nNeurons x kSequence matrix from network_spike_sequences
% PFpeaksSequence: sequence of Place field peaks
%
% Example usage, after running a simulation from randnet.m:
% [network_spike_sequences] = detect_PBE(spikes_V_m(network.E_indices,:), parameters);
% ranks_vec = network_spike_sequences(ithTrial).ranks_vec; % ranks for each detected PBE
% plot_PF_seq_corr(ranks_vec, PFpeaksSequence, network);


%% Default parameters:
seed = randi(10^6);
ithTrial = 1;
correlationType = 'Pearson'; % Pearson, Kendall, Spearman
usRelRank = 1;
nShuffles=1000;

%% Read in optional parameters, to overwrite above defaults
for i=1:2:length(varargin)
    switch varargin{i}
        case 'seed'
            seed = varargin{i+1};
        case 'ithTrial'
            ithTrial = varargin{i+1};
        case 'correlationType'
            correlationType = varargin{i+1};
        case 'usRelRank'
            usRelRank = varargin{i+1};
        case 'nShuffles'
            nShuffles = varargin{i+1};
        otherwise
            error('plot_clusterRates: Unknown input')
    end
end


%%

rng(seed)

if usRelRank
    ranks_vec = ranks_vec./sum(~isnan(ranks_vec), 1); % normalize to relative rank
end
y = PFpeaksSequence;

n = size(ranks_vec, 2);
r_actual=zeros(1, n);
p_actual=zeros(1, n);
for i=1:n
    [r_actual(i),p_actual(i)] = corr(ranks_vec(:,i),y,'type',correlationType, 'rows','complete');
end

r_shuff=zeros(1, nShuffles);
p_shuff=zeros(1, nShuffles);
for i=1:nShuffles
    x_shuff = ranks_vec( randperm(size(ranks_vec, 1)) , randi(size(ranks_vec, 2)) );
    %shuffID = randi(size(x, 2)); x_shuff = x( randperm(size(x, 1)) , shuffID)./sum(~isnan(x(:,shuffID)), 1);
    [r_shuff(i),p_shuff(i)] = corr(x_shuff,y,'type',correlationType, 'rows','complete');
end

figure; histogram(p_shuff, 10); title('Actual sequences')
xlabel('Correlation to PF (p-val)'); ylabel('Sequence (count)');

figure; histogram(p_actual, 10); title('Shuffled sequences')
xlabel('Correlation to PF (p-val)'); ylabel('Sequence (count)');

f = figure; hold on
histogram(r_shuff,  'BinWidth', 0.01, 'Normalization','probability', 'DisplayName','Shuffle')
histogram(r_actual, 'BinWidth', 0.01, 'Normalization','probability', 'DisplayName','Actual')
xlabel('Correlation to PF (r)')
ylabel('Sequences (probability)')
title([correlationType, ' correlations to PF sequence'])
legend('Location', 'Best')

end