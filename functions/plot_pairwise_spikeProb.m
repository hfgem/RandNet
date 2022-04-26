function plot_pairwise_spikeProb(spikes_V_m, PFpeaksSequence, varargin)
% Calculates and then plots data on the pairwaise spike probability--given
% a spike from one neuron, what is the probability that the subsequent spike
% will come from each of the other neurons
% 
% spikes_V_m: binary spike matrix from one trial
% PFpeaksSequence: sequence of place field peaks
%
% Example usage, after running a simulation from randnet.m:
% plot_pairwise_spikeProb(spikes_V_m(network.E_indices,:), PFpeaksSequence)
%

% TODO:
% Should only upper/lower triangle be considered?
% Should threshold rates be applied? 
% Should we only consider the single most likely to follow pairs?
%
% Note: bias towards upper triangle because there are many simultaneous
% spikes (see spikes_t) and rank of simultaneous spikes is determined by
% index in network structure
% Visualize simultaneous spikes with figure; histogram(spikes_t, max(spikes_t))



%% Default parameters:


%% Read in optional parameters, to overwrite above defaults
for i=1:2:length(varargin)
    switch varargin{i}
        otherwise
            error('plot_clusterRates: Unknown input')
    end
end


%%


[spikes_x,spikes_t] = find(spikes_V_m);


% For each cell, what is the probability that its spike is followed by a spike from each other cell
op = zeros(size(spikes_V_m, 1), size(spikes_V_m, 1)); 
for ithSpike = 1:(numel(spikes_t)-1)
    op(spikes_x(ithSpike), spikes_x(ithSpike+1) ) = op(spikes_x(ithSpike), spikes_x(ithSpike+1) ) + 1;
end


normOP = op./sum(op); normOP(isnan(normOP))=0;

figure; imagesc(normOP); colorbar
xlabel('next spike is from cell j'); ylabel('given spike from cell i'); title('Normalized probability')
figure; histogram(normOP); set(gca, 'YScale', 'log')
xlabel('Probability'); ylabel('Count (cell pair)'); title('Normalized probability')


figure; imagesc(op);
xlabel('next spike is from cell j'); ylabel('given spike from cell i'); title('Raw counts')
figure; histogram(op); set(gca, 'YScale', 'log')
xlabel('Post-synaptic cell spiked next (count)'); ylabel('Count (cell pair)'); title('Raw counts')


[B, I] = sort(op, 2, 'descend');
% figure; imagesc(B)
figure; errorbar(mean(B), var(B),'k-o', 'MarkerEdgeColor','r', 'MarkerSize', 3, 'CapSize',0); ylim([0, inf])
xlabel('Follower cell (sorted)'); ylabel('Spike pair count'); title('Mean (and var) over cells')

[B, Inorm] = sort(normOP, 2, 'descend');
% figure; imagesc(B)
figure; errorbar(mean(B), var(B),'k-o', 'MarkerEdgeColor','r', 'MarkerSize', 3, 'CapSize',0); ylim([0, inf])
xlabel('Follower cell (sorted)'); ylabel('Spike pair probability'); title('Mean (and var) over cells')


% Compare spike probability to PF rank order
x = PFpeaksSequence; PFpdist = squareform(pdist(x))./size(x, 1);
%PFpdist(isnan(PFpdist)) = 0;
figure; scatter( normOP(:), PFpdist(:))
xlabel('Probability of subsequent spike'); ylabel('Relative PF rank distance'); title('All cell pairs')
mdl = fitlm( normOP(:), PFpdist(:))


x = PFpeaksSequence; PFpdist = squareform(pdist(x))./size(x, 1);
%PFpdist(isnan(PFpdist)) = 0;
figure; scatter( op(:), PFpdist(:))
xlabel('Count of subsequent spike pair'); ylabel('Relative PF rank distance'); title('All cell pairs')
mdl = fitlm( op(:), PFpdist(:))


figure; scatter( Inorm(1,:), PFpeaksSequence)
xlabel('Most likely cell to follow'); ylabel('Cells place field rank'); title('Scatter of all E-cells')
mdl = fitlm(Inorm(1,:), PFpeaksSequence)

[Bpf, Ipf] = sort(PFpeaksSequence); Ipf(isnan(Bpf)) = nan;
figure; scatter( Inorm(1,:), Ipf)
xlabel('Most likely cell to follow'); ylabel('Sorted cells place field rank'); title('Scatter of all E-cells')
mdl = fitlm(Inorm(1,:), Ipf)


end