function [figHandleSubplots, figHandleSingle] = plot_dimRedNet(network, dimRedInput, varargin)
% Uses dimensionality reduction techniques to visualize the clustering of
% the network structure
%
% newtork: network structure
% dimRedInput: option for the input to the dimensionality reduction algorithms
%   'W': just the connection matrix
%   'WW': the connection matrix and its transpose [W, W']
%   'clust': for each neuron, the numer of connections it has to neurons in
%   each of the clustesr
%   'normClust': a neuron-wise normalization of clust
% Optional inputs:
%   seed: rng seed
%   E_only: if 1, only plot E-E connections
%   scatterSize: scatterSize(1) is the minimum size of the scatter points,
%   scattersize(2) is the scale for increasing the size of the scatter
%   points with cluster membership
%
% output: figure handles
%
% Example use, after running randnet.m:
% seed = 1; E_only = 1; dimRedInput = 'clust'; scatterSize = [1, 1];
% plot_dimRedNet(network, dimRedInput, 'seed', seed, 'E_only', E_only, 'scatterSize', scatterSize);

%% Defaults for optional parameters
seed =  randi(10^6); % rng seed, to replicate network construction
E_only = 1; % only plot E-E connections
scatterSize = [1, 1];

%% Read in optional parameters, to overwrite above defaults
for i=1:2:length(varargin)
    switch varargin{i}
        case 'seed'
            seed = varargin{i+1};
        case 'E_only'
            E_only = varargin{i+1};
        case 'scatterSize'
            scatterSize = varargin{i+1};
        otherwise
            error('plot_dimRedNet: Unknown input')
    end
end


%% Initialize seed and network, W, and plotting options
rng(seed)
if E_only
    network_indices = network.E_indices;
else
    network_indices = network.all_indices;
end

W = network.conns(network_indices, network_indices); % network conectivity matrix

c1 = network.cluster_mat(1,network_indices)';
c2 = network.cluster_mat(2,network_indices)';
c = c1*[0, 1, 0] + c2*[1, 0, 0]; % color of scatter plot points
sz = scatterSize(1) + scatterSize(2)*sum(network.cluster_mat(:,network_indices))'; % size of scatter plot points


%% set up matrix X, that will undergo dimensionality reduction
switch dimRedInput
    case 'W'
        X = [W]; % to use just outputs
    case 'WW'
        X = [W, W']; % to use both inputs and outputs
    case 'normClust'
        op = zeros(size(network.cluster_mat, 1), size(network.E_indices, 2));
        clustsE = network.cluster_mat(:,network.E_indices);
        for ithCell = 1:numel(network.E_indices)
            cellInd = network.E_indices(ithCell);
            inConns = network.conns(network.E_indices,cellInd);
            outConns = network.conns(cellInd,network.E_indices);
            inputClustCount =  clustsE*inConns;
            outputClustCount =  clustsE*outConns';
            op(:,ithCell) = inputClustCount + outputClustCount;
        end
        X = [op./sum(op)]';
    case 'clust'
        op = zeros(size(network.cluster_mat, 1), size(network.E_indices, 2));
        clustsE = network.cluster_mat(:,network.E_indices);
        for ithCell = 1:numel(network.E_indices)
            cellInd = network.E_indices(ithCell);
            inConns = network.conns(network.E_indices,cellInd);
            outConns = network.conns(cellInd,network.E_indices);
            inputClustCount =  clustsE*inConns;
            outputClustCount =  clustsE*outConns';
            op(:,ithCell) = inputClustCount + outputClustCount;
        end
        X = op';
    otherwise
        error('plot_dimRedNet: Unknown dimRedInput option')
end


%% Main function

Y_tsne = tsne(X);
try
    Y_mds = mdscale( pdist(X) ,2);
catch
    Y_mds = [nan nan];
end
[~,Y_PCA,~] = pca(X);


% Plot all three methods
figHandleSubplots = figure; 
sgtitle(['dim. red. plots of [' dimRedInput, ']'])
subplot(1,3,1);
plot(digraph(W), 'XData', Y_tsne(:,1),'YData', Y_tsne(:,2), ...
    'MarkerSize',sz, 'NodeColor', c, 'EdgeColor', 0.3*[1 1 1], 'EdgeAlpha',0.2)
title('tsne');

subplot(1,3,2); 
plot(digraph(W), 'XData', Y_mds(:,1),'YData', Y_mds(:,2), ...
    'MarkerSize',sz, 'NodeColor', c, 'EdgeColor', 0.3*[1 1 1], 'EdgeAlpha',0.2)
title('mds');

subplot(1,3,3);
plot(digraph(W), 'XData', Y_PCA(:,1),'YData', Y_PCA(:,2), ...
    'MarkerSize',sz, 'NodeColor', c, 'EdgeColor', 0.3*[1 1 1], 'EdgeAlpha',0.2)
title('PCA');


% Plot just t-sne
figHandleSingle = figure; 
Xpos = Y_tsne(:,1); Ypos = Y_tsne(:,2);
plot(digraph(W), 'XData', Xpos,'YData', Ypos, ...
    'MarkerSize',sz, 'NodeColor', c, 'EdgeColor', 0.3*[1 1 1], 'EdgeAlpha',0.2)
title(['tsne plot of [' dimRedInput, ']'])
hold on; h = zeros(4, 1);
h(1) = scatter(NaN,NaN,'r', 'filled');
h(2) = scatter(NaN,NaN,'g', 'filled');
h(3) = scatter(NaN,NaN,'y', 'filled');
h(4) = scatter(NaN,NaN,'k', 'filled');
legend(h, ' 1, ~2', '~1,  2', ' 1,  2', '~1, ~2', 'Location', 'Best');


% Plot along normalized dimensions
if isequal(dimRedInput, 'clust')
    % Based on normalized op, plot proximity to +1,-1 on plane 
    xVals = X(:,1)-X(:,2);
    yVals = X(:,3)-X(:,4);

    figure; plot(digraph(W), 'XData', xVals,'YData', yVals, ...
        'MarkerSize',sz, 'NodeColor', c, 'EdgeColor', 0.3*[1 1 1], 'EdgeAlpha',0.2)
    hold on; h = zeros(4, 1);
    h(1) = scatter(NaN,NaN,'r', 'filled');
    h(2) = scatter(NaN,NaN,'g', 'filled');
    h(3) = scatter(NaN,NaN,'y', 'filled');
    h(4) = scatter(NaN,NaN,'k', 'filled');
    legend(h, ' 1, ~2', '~1,  2', ' 1,  2', '~1, ~2', 'Location', 'Best');
    title('Norm. connection bias');
    xlabel('1-2 group connections'); ylabel('3-4 group connections'); 
end


end