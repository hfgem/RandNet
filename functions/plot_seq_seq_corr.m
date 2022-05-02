function plot_seq_seq_corr(ranks_vec, varargin)
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
nShuffles= min(2*size(ranks_vec, 2), 100);
maxNClust = round(sqrt( size(ranks_vec, 2) ));


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
        case 'nClust'
            maxNClust = varargin{i+1};
        otherwise
            error('plot_clusterRates: Unknown input')
    end
end

keyboard
%% Main:
rng(seed)

if usRelRank
    x_actual = ranks_vec./max(ranks_vec, [], 1);
else
    x_actual = ranks_vec;
end

actualRmat = corr(x_actual, 'type', correlationType, 'rows','pairwise');

% Plot sequence ranks
figure; imagesc(x_actual, 'AlphaData', ~isnan(x_actual)); colorbar
title('Sequence ranks'); xlabel('i^{th} sequence'); ylabel('Neuron');
set(gca,'Color','k')

% Plot clustered sequence ranks
square_pdist = squareform(pdist(actualRmat)); l = linkage(square_pdist); clustIDs = cluster(l, 'Maxclust',maxNClust); [~, I] = sort(clustIDs);
% figure; imagesc(x_actual(:,I)); colorbar
figure; imagesc(x_actual(:,I), 'AlphaData', ~isnan(x_actual(:,I))); colorbar
title('clustered rel. ranks'); xlabel('i^{th} sequence'); ylabel('Neuron');
set(gca,'Color','k')

% Plot clustered sorted (by one cluster) sequence ranks
[~, largestClustInd] = max(histcounts(clustIDs,'BinMethod','integers'));
meanRankLargestClust = nanmean(x_actual(:,clustIDs==largestClustInd), 2);
[a, cellRankCluster] = sort(meanRankLargestClust);
% figure; imagesc(x_actual(cellRankCluster,I)); colorbar
figure; imagesc(x_actual(cellRankCluster,I), 'AlphaData', ~isnan(x_actual(cellRankCluster,I))); colorbar
title('clustered, sorted rel. ranks'); xlabel('i^{th} sequence'); ylabel('Neuron');
set(gca,'Color','k')

% Plot sorted (neurons) sequence ranks
meanRankLargestClust = nanmean(x_actual, 2);
[a, cellRank] = sort(meanRankLargestClust);
% figure; imagesc(x_actual(cellRank,:)); colorbar
figure; imagesc(x_actual(cellRank,:), 'AlphaData', ~isnan(x_actual(cellRank,:))); colorbar
title('sorted (neurons) rel. ranks'); xlabel('i^{th} sequence'); ylabel('Neuron');
set(gca,'Color','k')

% Plot sorted (neurons) sorted (ranks) sequence ranks
meanRankLargestClust = nanmean(x_actual, 2);
[a, cellRank] = sort(meanRankLargestClust);
sequenceMeanCorr = zeros(size(x_actual, 2), 1);
for i = 1:numel(sequenceMeanCorr)
    x_temp = x_actual;
    %x_temp(x_actual==0) = nan;
    [~, eventSequence] = sort(x_temp(:,i));
    sequenceMeanCorr(i) = corr(eventSequence, cellRank, 'rows','complete');
end
[~, sequenceRank] = sort(sequenceMeanCorr);
% figure; imagesc(x_actual(cellRank,sequenceRank)); colorbar
figure; imagesc(x_actual(cellRank,sequenceRank), 'AlphaData', ~isnan(x_actual(cellRank,sequenceRank))); colorbar
title('sorted (neurons) sorted rel. ranks'); xlabel('i^{th} sequence'); ylabel('Neuron');
set(gca,'Color','k')


% perform dim. red. on relative ranks
ranks_vec(isnan(ranks_vec))=0;
Y_tsne = tsne(ranks_vec);
try
    Y_mds = mdscale( pdist(ranks_vec), 2); % this fails if any x_shuffle is identical to any x_actual
catch
    Y_mds = mdscale( pdist(ranks_vec), 2, 'start', 'random');
end
[~,Y_PCA,~] = pca(ranks_vec);

c = cellstr(num2str([1:size(ranks_vec, 1)]'));
c(3:end) = {''};
dx = 0.05; dy = 0.05; % displacement so the text does not overlay the data points
dx = 0.0; dy = 0.0; % displacement so the text does not overlay the data points

figure; sgtitle('relative rank dim. red.')
subplot(1,3,1); hold on; scatter(Y_tsne(:,1),Y_tsne(:,2)); text(Y_tsne(:,1)+dx, Y_tsne(:,2)+dy, c); title('tsne');
subplot(1,3,2); hold on; scatter(Y_mds(:,1),Y_mds(:,2)); text(Y_mds(:,1)+dx, Y_mds(:,2)+dy, c); title('mds');
subplot(1,3,3); hold on; scatter(Y_PCA(:,1),Y_PCA(:,2)); text(Y_PCA(:,1)+dx, Y_PCA(:,2)+dy, c); title('PCA');
% subplot(1,3,3); hold on; scatter(Y_PCA(:,2),Y_PCA(:,3)); text(Y_PCA(:,2)+dx, Y_PCA(:,3)+dy, c); title('PCA');


% Plot sorted be each cluster
figure
subplotCounter = 1;
nActualClusters = numel(unique(clustIDs));
for i = 1:nActualClusters
    
    meanRankLargestClust = nanmean(x_actual(:,clustIDs==i), 2);
    [a, cellRankCluster] = sort(meanRankLargestClust);
    
    % figure; imagesc(x_actual(cellRankCluster,I)); colorbar
    subplot(1,nActualClusters,subplotCounter); 
    imagesc(x_actual(cellRankCluster,I), 'AlphaData', ~isnan(x_actual(cellRankCluster,I))); 
    % colorbar
    %title('clustered, sorted rel. ranks'); xlabel('i^{th} sequence'); ylabel('Neuron');
    set(gca,'Color','k')

    subplotCounter = subplotCounter + 1;
end


% Plot ithClustxithClust sortings
figure
subplotCounter = 1;
nActualClusters = numel(unique(clustIDs));
for i = 1:nActualClusters
    
    meanRankLargestClust = nanmean(x_actual(:,clustIDs==i), 2);
    [a, cellRankCluster] = sort(meanRankLargestClust);
    
    for j = 1:nActualClusters
        
        % figure; imagesc(x_actual(cellRankCluster,I)); colorbar
        subplot(nActualClusters,nActualClusters,subplotCounter); 
        imagesc(x_actual(cellRankCluster,clustIDs==j), 'AlphaData', ~isnan(x_actual(cellRankCluster,clustIDs==j))); 
        % colorbar
        % title('clustered, sorted rel. ranks'); xlabel('i^{th} sequence'); ylabel('Neuron');
        set(gca,'Color','k')
        
        subplotCounter = subplotCounter + 1;
    end
end


%% Calculate inter-cluster correlations

nClusts = numel(unique(clustIDs));
op = zeros(numel(unique(clustIDs)), 1);
for i = 1:nClusts
    clustCorrs = actualRmat(clustIDs==i,clustIDs==i);
    clustCorrs(clustCorrs==1)=nan;
    op(i) = nanmean(clustCorrs, 'all');
end
histcounts(clustIDs,'BinMethod','integers')
intraClustCorrs = op'

temp = actualRmat;
temp(ismembertol(temp, 1, 10^-12))=nan; % remove self-correlations before nanmean
meanCorr = nanmean(temp, 'all')


%% Same as above, but for a shuffle
x_shuff = zeros( size(x_actual, 1), nShuffles);
for i = 1:nShuffles
    
    % randomly select from an actual sequence
    randSeq = x_actual(:, randi(size(x_actual, 2))); 
    firedInd = find(~isnan(randSeq));
    
    % Randomly permute only those cells that actually fired
    shufSeq = nan(size(randSeq));
    shufSeq(firedInd) = randSeq(firedInd(randperm(numel(firedInd))));
    
    x_shuff(:,i) = shufSeq;
end
shuffledRmat = corr(x_shuff, 'type', correlationType, 'rows','pairwise');

square_pdist = squareform(pdist(shuffledRmat)); l = linkage(square_pdist); clustIDs_shuff = cluster(l, 'Maxclust',maxNClust); [~, I] = sort(clustIDs_shuff);
nClusts = numel(unique(clustIDs_shuff));
shuffop = zeros(numel(unique(clustIDs_shuff)), 1);
for i = 1:nClusts
    clustCorrs = shuffledRmat(clustIDs_shuff==i,clustIDs_shuff==i);
    clustCorrs(ismembertol(clustCorrs, 1, 10^-12))=nan;
    shuffop(i) = nanmean(clustCorrs, 'all');
end
histcounts(clustIDs_shuff,'BinMethod','integers')
intraClustCorrs_shuff = shuffop'

temp = shuffledRmat;
temp(ismembertol(temp, 1, 10^-12))=nan;
overallMeanCorr_shuff = nanmean(temp, 'all')


% Plot clustered sorted (by one cluster) sequence ranks
[~, largestClustInd] = max(histcounts(clustIDs_shuff,'BinMethod','integers'));
meanRankLargestClust = nanmean(x_shuff(:,clustIDs_shuff==largestClustInd), 2);
[a, cellRankCluster] = sort(meanRankLargestClust);
% figure; imagesc(x_actual(cellRankCluster,I)); colorbar
figure; imagesc(x_shuff(cellRankCluster,I), 'AlphaData', ~isnan(x_shuff(cellRankCluster,I))); colorbar
title('clustered, sorted rel. ranks'); xlabel('i^{th} sequence'); ylabel('Neuron');
set(gca,'Color','k')


end