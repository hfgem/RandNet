function plot_meanRelRank(ranks_vec, PFpeaksSequence, varargin)
% calculate mean relative rank for actual and shuffle sequences, and then
% calcualtes the correlation of the mean sequence to the place field
% sequence
% 
% ranks_vec: the nENeurons x kSequence matrix from network_spike_sequences
% PFpeaksSequence: 1 x nENeurons sequence of place field peaks
%
% Example usage, after running a simulation from randnet.m:
% [network_spike_sequences] = detect_PBE(spikes_V_m(network.E_indices,:), parameters);
% ranks_vec = network_spike_sequences(ithTrial).ranks_vec;
% plot_meanRelRank(ranks_vec, PFpeaksSequence)


%% Default parameters:
seed = randi(10^6);
ithTrial = 1;
correlationType = 'Pearson'; % Pearson, Kendall, Spearman
minPartic = 1; % minimum number of sequences a cell needs to participate in


%% Read in optional parameters, to overwrite above defaults
for i=1:2:length(varargin)
    switch varargin{i}
        case 'seed'
            seed = varargin{i+1};
        case 'ithTrial'
            ithTrial = varargin{i+1};
        case 'correlationType'
            correlationType = varargin{i+1};
        case 'minPartic'
            minPartic = varargin{i+1};
        otherwise
            error('plot_clusterRates: Unknown input')
    end
end


%% Main
rng(seed)

relPFRank = PFpeaksSequence./size(PFpeaksSequence, 1) ;  
ranks_vec = ranks_vec./sum(~isnan(ranks_vec), 1); % normalize to relative rank

meanRelRank = nanmean(ranks_vec, 2);
stdRelRank = nanstd(ranks_vec, [], 2);
meanRelRank(sum(~isnan(ranks_vec), 2)<=minPartic) = nan;

[r,p_actual] = corr(meanRelRank,relPFRank,'type',correlationType, 'rows','complete')
[B,I_actual] = sort(meanRelRank);

figure; hold on; 
errorbar(meanRelRank(I_actual), 1:numel(meanRelRank),stdRelRank(I_actual)./sqrt(size(ranks_vec, 2)),'horizontal', 'CapSize',0)
plot(meanRelRank(I_actual), 1:numel(meanRelRank), 'k')
xlabel('Relative rank (mean \pm SEM)'); ylabel('Neuron (sorted)'); 
title(['Actual sequences: p=', num2str(p_actual)]);
hold on; scatter(relPFRank(I_actual), 1:numel(relPFRank))

% Same mean rank as above, but for x_shuff
nShuf = size(ranks_vec, 2);
x_shuff = zeros( size(ranks_vec, 1), nShuf);
for i = 1:nShuf
    ithSeq = ranks_vec(:,i);
    firedInd = find( ~isnan(ithSeq) );
    
    shufSeq = nan(size(ithSeq));
    shufSeq(firedInd) = ithSeq(firedInd(randperm(numel(firedInd))));
    
    x_shuff(:,i) = shufSeq;
end

meanRelRank_shuff = nanmean(x_shuff, 2);
stdRelRank = nanstd(x_shuff, [], 2);
meanRelRank_shuff(sum(~isnan(x_shuff), 2)<=minPartic) = nan;

[r,p_shuffle] = corr(meanRelRank_shuff,relPFRank,'type',correlationType, 'rows','complete');
[B,I_shuffle] = sort(meanRelRank_shuff);

figure; hold on
errorbar(meanRelRank_shuff(I_shuffle), 1:numel(meanRelRank_shuff),stdRelRank(I_shuffle)./sqrt(size(x_shuff, 2)),'horizontal', 'CapSize',0)
plot(meanRelRank_shuff(I_shuffle), 1:numel(meanRelRank_shuff), 'k')
xlabel('Relative rank (mean \pm SEM)'); ylabel('Neuron (sorted)'); 
title(['Shuffled sequences: p=', num2str(p_shuffle)]);
hold on; scatter(relPFRank(I_actual), 1:numel(relPFRank))

end