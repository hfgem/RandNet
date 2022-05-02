function [f1Handle, f2Handle, f3Handle] = plot_PF_seq_corr(ranks_vec, PFpeaksSequence, varargin)
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


%% Main:

rng(seed)

if usRelRank
    ranks_vec = ranks_vec./sum(~isnan(ranks_vec), 1); % normalize to relative rank
end
y = PFpeaksSequence;

[r_actual,p_actual] = corr(ranks_vec, y, 'type', correlationType, 'rows', 'pairwise');

x_shuff = zeros( size(ranks_vec, 1), nShuffles); 
for i = 1:nShuffles
    % randomly select from an actual sequence
    randSeq = ranks_vec(:, randi(size(ranks_vec, 2))); 
    firedInd = find(~isnan(randSeq));
    
    % Randomly permute only those cells that actually fired
    shufSeq = nan(size(randSeq));
    shufSeq(firedInd) = randSeq(firedInd(randperm(numel(firedInd))));
    x_shuff(:,i) = shufSeq;
end

[r_shuff,p_shuff] = corr(x_shuff,y,'type',correlationType, 'rows', 'pairwise');


f2Handle = figure; 
histogram(p_shuff, 10); title('Shuffled sequences')
xlabel('Correlation to PF (p-val)'); ylabel('Sequence (count)');

f3Handle = figure; 
histogram(p_actual, 10); title('Actual sequences')
xlabel('Correlation to PF (p-val)'); ylabel('Sequence (count)');

f1Handle = figure; hold on
histogram(r_shuff,  'BinWidth', 0.01, 'Normalization','probability', 'DisplayName','Shuffle')
histogram(r_actual, 'BinWidth', 0.01, 'Normalization','probability', 'DisplayName','Actual')
xlabel('Correlation to PF (r)')
ylabel('Sequences (probability)')
title([correlationType, ' correlations to PF sequence'])
legend('Location', 'Best')

end