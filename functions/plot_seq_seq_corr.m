function plot_seq_seq_corr(ranks_vec, varargin)
% Calculates sequence to sequences correlations of the ranks_vec produced by 
% detect_PBE.
% For both actual and shuffled data, sequences are clustered and plotted.
% 
% ranks_vec: the nNeurons x kSequence matrix from network_spike_sequences
%
% Example usage, after running a simulation from randnet.m:
% [network_spike_sequences] = detect_PBE(E_spikes_V_m, parameters);
% ranks_vec = network_spike_sequences.ranks_vec;
% plot_seq_seq_corr(ranks_vec, 'correlationType', correlationType);


%% Default parameters:
seed = randi(10^6);
correlationType = 'Pearson'; % Pearson, Kendall, Spearman
usRelRank = 1;
nShuffles= max(2*size(ranks_vec, 2), 100);
maxNClust = round(sqrt( size(ranks_vec, 2) ));


%% Read in optional parameters, to overwrite above defaults
for i=1:2:length(varargin)
    switch varargin{i}
        case 'seed'
            seed = varargin{i+1};
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


%% Main:
rng(seed)

% Actual ranks and correlations
if usRelRank
    x_actual = ranks_vec./max(ranks_vec, [], 1);
else
    x_actual = ranks_vec;
end
actualRmat = corr(x_actual, 'type', correlationType, 'rows','pairwise');

% Shuffle ranks and correlations
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

% Combined cross-correlation matrix
combinedRmat = corr([x_actual, x_shuff], 'type', correlationType, 'rows','pairwise');

%% Clustering and Plotting

% Actual data: 
% Create clusters
square_pdist = squareform(pdist(actualRmat)); l = linkage(square_pdist); clustIDs = cluster(l, 'Maxclust',maxNClust); [~, clustSort_Actual] = sort(clustIDs);

% Plot relative ranks, neurons sorted by mean rank, sequences clusterd
meanRelRank_actual = nanmean(x_actual, 2);
[a, cellRank_actual] = sort(meanRelRank_actual);
% figure; imagesc(x_actual(cellRank,:)); colorbar
figure; imagesc(x_actual(cellRank_actual,clustSort_Actual), 'AlphaData', ~isnan(x_actual(cellRank_actual,clustSort_Actual))); colorbar
title('Clustered, sorted rel. ranks'); xlabel('i^{th} sequence'); ylabel('Neuron');
set(gca,'Color','k')

% Plot sorted be each cluster
figure
subplotCounter = 1;
nActualClusters = numel(unique(clustIDs));
for i = 1:nActualClusters
    meanRelRank_actual = nanmean(x_actual(:,clustIDs==i), 2);
    [a, cellRankCluster] = sort(meanRelRank_actual);
    
    subplot(1,nActualClusters,subplotCounter); 
    imagesc(x_actual(cellRankCluster,clustSort_Actual), 'AlphaData', ~isnan(x_actual(cellRankCluster,clustSort_Actual))); 
    set(gca,'Color','k')

    subplotCounter = subplotCounter + 1;
end
sgtitle('Clustered sequences, neurons sorted by each cluster')


% Shuffled data: 
% Create clusters
square_pdist = squareform(pdist(shuffledRmat)); l = linkage(square_pdist); clustIDs_shuff = cluster(l, 'Maxclust',maxNClust); [~, clustSort_shuff] = sort(clustIDs_shuff);

% Plot relative ranks, neurons sorted by mean rank, sequences clusterd
meanRelRank_shuffle = nanmean(x_shuff, 2);
[a, cellRank_shuffle] = sort(meanRelRank_shuffle);
% figure; imagesc(x_actual(cellRank,:)); colorbar
figure; imagesc(x_shuff(cellRank_shuffle,clustSort_shuff), 'AlphaData', ~isnan(x_shuff(cellRank_shuffle,clustSort_shuff))); colorbar
title('Shuf: Clustered, sorted rel. ranks'); xlabel('i^{th} sequence'); ylabel('Neuron');
set(gca,'Color','k')

% Plot sorted be each cluster
figure
subplotCounter = 1;
nActualClusters = numel(unique(clustIDs_shuff));
for i = 1:nActualClusters
    meanRelRank_shuffle = nanmean(x_shuff(:,clustIDs_shuff==i), 2);
    [a, cellRankCluster] = sort(meanRelRank_shuffle);
    
    subplot(1,nActualClusters,subplotCounter); 
    imagesc(x_shuff(cellRankCluster,clustSort_shuff), 'AlphaData', ~isnan(x_shuff(cellRankCluster,clustSort_shuff))); 
    set(gca,'Color','k')

    subplotCounter = subplotCounter + 1;
end
sgtitle('Shuf: Clustered rel. ranks, cluster-sorted')


x = actualRmat; p = squareform(pdist(x)); l = linkage(p); c = cluster(l, 'Maxclust',maxNClust); [~, I] = sort(c);
figure; imagesc(x(I,I)); colorbar
title('Clustered ixj Pearson Corr.'); xlabel('i^{th} sequence'); ylabel('j^{th} sequence');

x = shuffledRmat; p = squareform(pdist(x)); l = linkage(p); c = cluster(l, 'Maxclust',maxNClust); [~, I] = sort(c);
figure; imagesc(x(I,I)); colorbar
title('Shuf: Clustered ixj Pearson Corr.'); xlabel('i^{th} sequence'); ylabel('j^{th} sequence');


%% Plot Dimensionality reduction:

% Relative ranks:
dimRedData = [x_actual, x_shuff]';
dimRedData(isnan(dimRedData))=0;
Y_tsne = tsne(dimRedData);
try
    Y_mds = mdscale( pdist(dimRedData), 2); % this fails if any x_shuffle is identical to any x_actual
catch
    Y_mds = mdscale( pdist(dimRedData), 2, 'start', 'random');
end
[~,Y_PCA,~] = pca(dimRedData);

clr = [repmat([0 0 0], size(x_actual, 2), 1); repmat([1 0 0], size(x_shuff, 2), 1) ];
figure; sgtitle('relative rank dim. red. (red=shuf.)')
subplot(1,3,1); hold on; scatter(Y_tsne(:,1),Y_tsne(:,2), [], clr); title('t-sne');
subplot(1,3,2); hold on; scatter(Y_mds(:,1),Y_mds(:,2), [], clr); title('MDS');
subplot(1,3,3); hold on; scatter(Y_PCA(:,1),Y_PCA(:,2), [], clr); title('PCA');


% Sequence correlations:
dimRedData = [combinedRmat];
Y_tsne = tsne(dimRedData);
try
    Y_mds = mdscale( pdist(dimRedData), 2); % this fails if any x_shuffle is identical to any x_actual
catch
    Y_mds = mdscale( pdist(dimRedData), 2, 'start', 'random');
end
[~,Y_PCA,~] = pca(dimRedData);

clr = [repmat([0 0 0], size(x_actual, 2), 1); repmat([1 0 0], size(x_shuff, 2), 1) ]  ; size(clr)
figure; sgtitle('Sequence correlations dim. red. (red=shuf.)')
subplot(1,3,1); hold on; scatter(Y_tsne(:,1),Y_tsne(:,2), [], clr); title('t-sne');
subplot(1,3,2); hold on; scatter(Y_mds(:,1),Y_mds(:,2), [], clr); title('MDS');
subplot(1,3,3); hold on; scatter(Y_PCA(:,1),Y_PCA(:,2), [], clr); title('PCA');


%% Mean overal and inter-cluster correlations

% Actual:
nClusts = numel(unique(clustIDs));
op = zeros(numel(unique(clustIDs)), 1);
for i = 1:nClusts
    clustCorrs = actualRmat(clustIDs==i,clustIDs==i);
    clustCorrs(clustCorrs==1)=nan;
    op(i) = nanmean(clustCorrs, 'all');
end
clusterCounts_actual = histcounts(clustIDs,'BinMethod','integers')
intraClustCorrs = op'

temp = actualRmat;
temp(ismembertol(temp, 1, 10^-12))=nan; % remove self-correlations before nanmean
meanCorr = nanmean(temp, 'all')


% Shuffle:
nClusts = numel(unique(clustIDs_shuff));
shuffop = zeros(numel(unique(clustIDs_shuff)), 1);
for i = 1:nClusts
    clustCorrs = shuffledRmat(clustIDs_shuff==i,clustIDs_shuff==i);
    clustCorrs(ismembertol(clustCorrs, 1, 10^-12))=nan;
    shuffop(i) = nanmean(clustCorrs, 'all');
end
clusterCounts_shuffle = histcounts(clustIDs_shuff,'BinMethod','integers')
intraClustCorrs_shuff = shuffop'

temp = shuffledRmat;
temp(ismembertol(temp, 1, 10^-12))=nan;
overallMeanCorr_shuff = nanmean(temp, 'all')


end