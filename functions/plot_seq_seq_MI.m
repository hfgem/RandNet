function plot_seq_seq_MI(ranks_vec, varargin)
% Calculates sequence to sequences Matching Index of the ranks_vec produced by 
% detect_PBE.
% For both actual and shuffled data, sequences are clustered and plotted.
% 
% ranks_vec: the nNeurons x kSequence matrix from network_spike_sequences
%
% Example usage, after running a simulation from randnet.m:
% [network_spike_sequences] = detect_PBE(E_spikes_V_m, parameters);
% ranks_vec = network_spike_sequences(ithTrial).ranks_vec;
% plot_seq_seq_MI(ranks_vec)


%% Default parameters:
seed = randi(10^6);
nShuffles= max(2*size(ranks_vec, 2), 100);
maxNClust = round(sqrt( size(ranks_vec, 2) ));
penalize_nonspike = 0;


%% Read in optional parameters, to overwrite above defaults
for i=1:2:length(varargin)
    switch varargin{i}
        case 'seed'
            seed = varargin{i+1};
        case 'nShuffles'
            nShuffles = varargin{i+1};
        case 'maxNClust'
            maxNClust = varargin{i+1};
        case 'penalize_nonspike'
            penalize_nonspike = varargin{i+1};
        otherwise
            error('plot_clusterRates: Unknown input')
    end
end


%% Main
rng(seed)

X = ranks_vec; % ranks for each detected sequence
X_n = isnan(ranks_vec); % non-spiking neurons, for each sequence

% MI, actual
[matching_index, matching_index_mod] = calculate_trajectory_similarity_mi2(X, ...
    X_n, penalize_nonspike);

x = matching_index; p = squareform(pdist(x)); l = linkage(p); c = cluster(l, 'Maxclust',maxNClust); [~, I] = sort(c);
figure; imagesc(x(I,I)); colorbar
title('Clustered ixj MI'); xlabel('i^{th} sequence'); ylabel('j^{th} sequence');


% MI, shuffle
x_shuff = zeros( size(X, 1), nShuffles);
for i = 1:nShuffles
    % randomly select from an actual sequence
    randSeq = X(:, randi(size(X, 2))); 
    firedInd = find(~isnan(randSeq));
    % Randomly permute only those cells that actually fired
    shufSeq = nan(size(randSeq));
    shufSeq(firedInd) = randSeq(firedInd(randperm(numel(firedInd))));
    x_shuff(:,i) = shufSeq;
end
x_shuff_n = isnan(x_shuff);

[matching_indexShuff, matching_index_modShuff] = calculate_trajectory_similarity_mi2(x_shuff, ...
    x_shuff_n, penalize_nonspike);

x = matching_indexShuff; p = squareform(pdist(x)); l = linkage(p); c = cluster(l, 'Maxclust',maxNClust); [~, I] = sort(c);
figure; imagesc(x(I,I)); colorbar
title('Shuffle: Clustered ixj MI'); xlabel('i^{th} sequence'); ylabel('j^{th} sequence');


% Histogram:
nantri = tril(ones(size(matching_index)), -1); nantri(nantri==0) = nan;
nantriShuff = tril(ones(size(matching_indexShuff)), -1); nantriShuff(nantriShuff==0) = nan;

f = figure; hold on
histogram(matching_indexShuff.*nantriShuff, 'BinWidth', 0.01, 'Normalization','probability', 'DisplayName','Shuffle')
histogram(matching_index.*nantri, 'BinWidth', 0.01, 'Normalization','probability', 'DisplayName','Actual')
xlabel('Matching index'); 
ylabel('sequence pair (count)');
title('MI between all sequence pairs')
legend('Location', 'Best')


end

