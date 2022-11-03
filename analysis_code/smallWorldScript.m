% smallWorldscript.m
% Initial script section from PFoptSeq_scriptv1_0.m

%% Calculate example network's small-world statistic
%
% C: clustering coefficient
% L: mean path length
%
% n: network size (?)
% k: mean number of connections for each node (?)
% 
% _r: random graph
% _l: lattice graph

% N = 1000; K = 10; beta = 0.005;
% N = numel(network.E_indices); beta = 0.07; K = parameters.conn_prob*numel(network.E_indices); 

% W = network.conns(network.E_indices,network.E_indices); % E-E connection matrix
% W = rand(size(network.conns)) < parameters.conn_prob;
% W = adjacency(WattsStrogatz(N,K,beta));

N = 500;
conProb = 0.08;
W = [rand(N)<conProb];

n = size(W, 1)
k = sum(W, 'all')/n % divide by 2 if symmetric?


% Note: not sure of _r and _l are correct for digraphs
% https://www.cambridge.org/core/services/aop-cambridge-core/content/view/56ED3E109EE91CEA6FAE2D82FF20DB30/S2050124217000054a.pdf/how-small-is-it-comparing-indices-of-small-worldliness.pdf
C_r = k/n
L_r = log(n)/log(k)

C_l = (3 * (k-2) )/(4 * (k-1) )
L_l = n / (2*k)

L = mean(distances(digraph(W)), 'all')

[C, ~] = avgClusteringCoefficient(W)

% Small world index, ->1 indicates small-world like
SWI_actual = [(L - L_l) / (L_r - L_l)] * [(C - C_r) / (C_l - C_r) ]

% disp(parameters.clusters)
% disp(parameters.mnc)

%{
% Clustering coefficient for directed graph, eq 8 from
% https://journals.aps.org/pre/pdf/10.1103/PhysRevE.76.026107
A = W;
numer = (A+A')^3;
dtot = sum(A, 1) + sum(A', 1);
dbi = A^2;
for i = 1:size(W, 1)
    Cd_i(i) = [numer(i,i)] / [2*( dtot(i)*(dtot(i)-1) - (2*dbi(i,i))  )] ;
end
Cd = mean(Cd_i);
%}

%% Verify analytic solutions for C_r and L_r for directed graphs

nNets = 50;

N = 1000;
p = 0.03;

LVec = zeros(1, nNets);
CVec = zeros(1, nNets);

tic
for ithNet = 1:nNets
    
    W = [rand(N)<=p];
    % W = [rand(N)<=(p/2)]; W = W+W';

    n = size(W, 1);
    k = sum(W, 'all')/n; % divide by 2 if symmetric?

    wDistances = distances(digraph(W));
    wDistances(wDistances==0)=nan;
    L = nanmean(wDistances, 'all');
    
    C = fagioloCC(W);
    %C= avgClusteringCoefficient(W);
    
    LVec(ithNet) = L;
    CVec(ithNet) = C;
    
end
runTime = toc

kPred = N*p;
egamma = double(eulergamma);

L_r = ( (log(N)-egamma) /log(kPred)) + 0.5 % https://doi.org/10.48550/arXiv.cond-mat/0407098
mean(LVec)

C_r = p
mean(CVec)

figure; hold on; xline(L_r, 'r', 'LineWidth', 2); xline(mean(LVec), 'k', 'LineWidth', 2); histogram(LVec); 
title('L_r for random networks'); legend({'Analytic', 'Mean', 'Monte Carlo'}, 'Location', 'Best')

figure; hold on; xline(C_r, 'r', 'LineWidth', 2); xline(mean(CVec), 'k', 'LineWidth', 2); histogram(CVec);
title('C_r for random networks'); legend({'Analytic', 'Mean', 'Monte Carlo'}, 'Location', 'Best')


%% Verify analytic solutions for C_l and L_l for directed graphs

nNets = 1;

N = 400;
p = 0.125
beta = 0.00;

LVec = zeros(1, nNets);
CVec = zeros(1, nNets);

ithNet = 1;
W = adjacency(WattsStrogatz(N, N*p, beta));
% figure; imagesc(W); keyboard

n = size(W, 1);
k = sum(W, 'all')/n; % divide by 2 if symmetric?

wDistances = distances(digraph(W));
wDistances(wDistances==0)=nan;
L = nanmean(wDistances, 'all');

C = fagioloCC(W);
% [C, ~] = avgClusteringCoefficient(W);

LVec(ithNet) = L;
CVec(ithNet) = C;

egamma = double(eulergamma);

L_l = n / (2*k)
L_l = (n-egamma) / (2*k)+0.5
L_l = (n*(n+k-2) ) / (2*k*(n-1)) 
mean(LVec)

% C_r = kPred/N;
% C_r = 4*(N-1)*(N-2)* ((p/2)^3) % https://journals.aps.org/pre/pdf/10.1103/PhysRevE.76.026107
C_l = (3 * (k-2) )/(4 * (k-1) )
mean(CVec)


figure; hold on; xline(L_l, '--k', 'LineWidth', 2); xline(LVec, ':r', 'LineWidth', 2); 
title('L_l for lattice networks'); legend({'Analytic', 'Monte Carlo'}, 'Location', 'Best')

figure; hold on; xline(C_l, '--k', 'LineWidth', 2); xline(CVec, ':r', 'LineWidth', 2); 
title('C_l for lattice networks'); legend({'Analytic', 'Monte Carlo'}, 'Location', 'Best')



%% Calculate SWI for WS network
% WS network, an undirected graph with
% N: number of nodes
% K: mean degree
% beta: probability of randomly re-wiring an edge from an initial ring
%     lattice network
%
% (N*K)/2 number of edges

% Default parameters
 params.N = 1000; params.K = 200; params.beta = 0.01;
 params.N = 375; params.K = round(params.N* (0.08/2)); params.beta = 0.01;

% Parameters to vary
 variedParam.name = 'beta'; variedParam.range = 0:0.01:0.2; % for beta
% variedParam.name = 'K'; variedParam.range = 1:25; % for K
% variedParam.name = 'N'; variedParam.range = 100:100:2000; % for N

egamma = double(eulergamma);

tic
op = zeros(size(variedParam.range));
for ithParam = 1:numel(variedParam.range)
    
	params.(variedParam.name) = variedParam.range(ithParam);
    
    W = adjacency(WattsStrogatz(params.N,params.K,params.beta));

    n = size(W, 1);
    k = sum(W, 'all')/n; % divide by 2 if symmetric?

    % Note: not sure of _r and _l are correct for digraphs
    % https://www.cambridge.org/core/services/aop-cambridge-core/content/view/56ED3E109EE91CEA6FAE2D82FF20DB30/S2050124217000054a.pdf/how-small-is-it-comparing-indices-of-small-worldliness.pdf
    C_r = k/n;
    L_r = log(n)/log(k);
    C_l = (3 * (k-2) )/(4 * (k-1) );
    L_l = n / (2*k);
    
    % Neal, 2017. Note: not sure of _r and _l are correct for digraphs  
    %{
    k = sum(W, 'all')/n; % divide by 2 if symmetric?
    C_r = k/n;
    L_r = log(n)/log(k);
    C_l = (3 * (k-2) )/(4 * (k-1) );
    L_l = n / (2*k);
     %}

    % Validated calculations
    k = sum(W, 'all')/n;
    p = sum(W, 'all')/(n^2-n);
    C_r = p;
    L_r = ( (log(n)-egamma) /log(k)) + 0.5;
    C_l = (3 * (k-2) )/(4 * (k-1) );
    L_l = (n-egamma) / (2*k)+0.5;

    L = mean(distances(digraph(W)), 'all');

    [C, ~] = avgClusteringCoefficient(W);

    % Small world index, ->1 indicates small-world like
    SWI_actual = [(L - L_l) / (L_r - L_l)] * [(C - C_r) / (C_l - C_r) ];

    op(ithParam) = SWI_actual;
    
    %{
    figure; imagesc(W)
    %}
end
runTime = toc

otherParams = params; otherParams = rmfield(otherParams, variedParam.name);
figure; plot(variedParam.range, op)
xlabel([variedParam.name]); ylabel('SWI')
% title(strtrim(formattedDisplayText(otherParams)))
title(erase(formattedDisplayText(otherParams), newline))
ylim([0, 1])


%% Calculate SWI for ER random graphs

nNets = 100;

N = 375;
p = 0.08;

op = zeros(1, nNets);

tic
for ithNet = 1:nNets
    
    W = [rand(N)<=p]; % directed
    % W = [rand(N)<=(p/2)]; W = W+W'; % un-directed
    
    wDistances = distances(digraph(W));
    wDistances(wDistances==0)=nan;
    L = nanmean(wDistances, 'all');
    
    C = fagioloCC(W);
    %C= avgClusteringCoefficient(W);
    
    % Validated calculations
    n = size(W, 1);
    k = sum(W, 'all')/n;
    p = sum(W, 'all')/(n^2-n);
    C_r = p;
    L_r = ( (log(n)-egamma) /log(k)) + 0.5;
    C_l = (3 * (k-2) )/(4 * (k-1) );
    L_l = (n-egamma) / (2*k)+0.5;
    

    % Small world index, ->1 indicates small-world like
    SWI_actual = [(L - L_l) / (L_r - L_l)] * [(C - C_r) / (C_l - C_r) ];

    % if SWI_actual > 1; keyboard; end

    op(ithNet) = SWI_actual;
    
end
runTime = toc

figure; histogram(op); title('SWI, ER random graph')
% xlim([0, 0.1])


%% Calculate SWI for randNet across nclusters and mnc

rng(1)

addpath('~/Documents/GitHub/RandNet/functions')

parameters.conn_prob = 0.08;
parameters.n = 500;
%  mncVec = 1:0.5:21; nClustersVec = 1:1:21;
mncVec = 1:0.5:6;  nClustersVec = 1:1:36;
mncVec = 1:0.5:11; nClustersVec = 2:1:21; % to match 'results_2022-05-08T15-43.mat'
% mncVec = 1:0.25:6; nClustersVec = 2:2:42; % to match 'results_2022-04-29T23-05.mat'

mncVec = 1:0.25:3.0; nClustersVec = 5:5:60; % to match 'results_....mat'
mncVec = 1:0.025:3.0; nClustersVec = 1:1:60; % finer grid
% mncVec = 1.5:0.5:5.0; nClustersVec = 2:3:25; %


parameters.p_E = 0.75;
parameters.p_I = 1-parameters.p_E;
parameters.n_I = round((1-parameters.p_E)*parameters.n); %number of inhibitory neurons
parameters.n_E = parameters.n - parameters.n_I;
parameters.W_gin = 1; parameters.cueSigma = 1; parameters.IcueScale = 1;

parameters.envIDs = 1;
include_all = 3; % =3 means all neurons in at least one cluster

egamma = double(eulergamma);

% mncVec = 1:0.5:1; nClustersVec = 1:1:4;

tic
op = zeros( numel(mncVec), numel(nClustersVec) );
for ithmnc = 1:numel(mncVec)
    
    parameters.mnc = mncVec(ithmnc);
    
    for ithnClust = 1:numel(nClustersVec)
        
        parameters.clusters = nClustersVec(ithnClust);

        if parameters.mnc>parameters.clusters
            op(ithmnc, ithnClust) = nan;
            continue
        end
        
        parameters.cluster_n = round((parameters.mnc*parameters.n) / parameters.clusters) ; % Rather than above method, explicitly declare mnc as a parameter
        parameters.npairs = parameters.n*(parameters.n-1); %total number of possible neuron connections
        parameters.nclusterpairs = parameters.cluster_n*(parameters.cluster_n - 1)*parameters.clusters; %total number of possible intra-cluster connections
        parameters.cluster_prob = min(parameters.conn_prob*parameters.npairs/parameters.nclusterpairs,1); %0.2041; %intra-cluster connection probability
        
        network = create_clusters(parameters, 'include_all', include_all);
        W = network.conns(network.E_indices, network.E_indices);
        W = W>0;
        n = size(W, 1);        
        
        % Enforce symmetry in a randnet network
        %{
        parameters.conn_prob = 0.04;
        network = create_clusters(parameters);
        W = network.conns(network.E_indices, network.E_indices);
        W = W+W'; W = W>0;
        n = size(W, 1);
        %}
        
        % Neal, 2017. Note: not sure of _r and _l are correct for digraphs  
        %{
        k = sum(W, 'all')/n/2; % divide by 2 if symmetric?
        C_r = k/n;
        L_r = log(n)/log(k);
        C_l = (3 * (k-2) )/(4 * (k-1) );
        L_l = n / (2*k);
         %}
        
        % Validated calculations
        k = sum(W, 'all')/n;
        p = sum(W, 'all')/(n^2-n);
        C_r = p;
        L_r = ( (log(n)-egamma) /log(k)) + 0.5;
        C_l = (3 * (k-2) )/(4 * (k-1) );
        L_l = (n-egamma) / (2*k)+0.5;
       
        %{
        WLattice = adjacency(WattsStrogatz(n,round((p/2)*n),0));
        % C_l = fagioloCC(WLattice);
        wDistances_lattice = distances(digraph(WLattice));
        wDistances_lattice(wDistances_lattice==0)=nan;
        L_l = nanmean(wDiwDistances_latticestances, 'all');        
        %}
        
        wDistances = distances(digraph(W));
        wDistances(wDistances==0)=nan;
        L = nanmean(wDistances, 'all');
    
        %[C, ~] = avgClusteringCoefficient(W);
        C = fagioloCC(W);
        % ci = clusteringcoef(W); C = mean(ci);
        
        % Small world index, ->1 indicates small-world like
        SWI_actual = [(L - L_l) / (L_r - L_l)] * [(C - C_r) / (C_l - C_r) ];
        
         % if SWI_actual > 1; disp( [(L - L_l) / (L_r - L_l)]); disp([(C - C_r) / (C_l - C_r) ]); keyboard; end
        
        op(ithmnc, ithnClust) = SWI_actual;
        
        %{
        
        % figure; imagesc(W) % Plot connectivity matrix
        maxNClust = parameters.clusters;
        x = W; p = squareform(pdist(x)); l = linkage(p); c = cluster(l, 'Maxclust',maxNClust); [~, I] = sort(c);
        figure; imagesc(x(I,I)); colorbar
        title(['Clustered W, nClust=', num2str(maxNClust)]); xlabel('Neuron'); ylabel('Neuron');
        
        %}
        
        %{
        % Debuging lines
        figure; histogram(sum(network.cluster_mat(:,network.E_indices), 1))
        parameters.mnc
        parameters.clusters
        SWI_actual
        keyboard
        %}
        

    end
end
runTime = toc

wrongConnectivity = parameters.n_E*[(1./nClustersVec).*mncVec'] < parameters.n_E*parameters.conn_prob;
wrongConnectivity = ([(1./nClustersVec).*(parameters.n_E.*mncVec')].*mncVec') < parameters.n_E*parameters.conn_prob;

op2 = op;
op2(wrongConnectivity) = nan;
AlphaData =  (~isnan(op2')+0.25) / 1.5;

figure; imagesc(mncVec, nClustersVec, op', 'AlphaData', AlphaData);
colorbar; set(gca,'YDir','normal')
title('SWI'); xlabel('mnc'); ylabel('nClusters')

% halfPerctlEdge = 20; caxis([prctile(op, halfPerctlEdge, 'all'), prctile(op, 100-halfPerctlEdge, 'all')])

% caxis([0, inf])
% caxis([prctile(op, 40, 'all'), prctile(op, 95, 'all')])
% caxis([prctile(op, 20, 'all'), prctile(op, 50, 'all')])


% Detect threshold crossing, then plot fitted exponential curve
plotThreshold = 1;
if plotThreshold
    thresh = 0.4;
    x = op2>thresh; 
    % figure; imagesc(mncVec, nClustersVec, x'); set(gca, 'YDir', 'normal')
    % figure; imagesc(mncVec, nClustersVec,  diff(x, [], 2)'); set(gca, 'YDir', 'normal')
    [xvalinds, yvalinds] = find( diff(x, [], 2)==1 ) ;

    ft = fittype('a*exp(b*x) + c');
    % [fexp_offset, fexp_offset_G] = fit(mncVec(xvalinds)', nClustersVec(yvalinds)', ft);
    [fexp, fexp_G] = fit(mncVec(xvalinds)', nClustersVec(yvalinds)','exp1');

    hold on; scatter(mncVec(xvalinds), nClustersVec(yvalinds), 'r')
    plot(mncVec, fexp(mncVec), 'r')
end

% Additional plots related to above
plotExtra = 1; 
if plotExtra
    
    % How many neurons are in each cluster
    nNeuronsPerCluster = [(1./nClustersVec).*(parameters.n_E.*mncVec')];
    figure; imagesc(mncVec, nClustersVec, nNeuronsPerCluster'); %, 'AlphaData', ~isnan(op'));
    colorbar; set(gca,'YDir','normal')
    title('nNeuronsPerCluster'); xlabel('mnc'); ylabel('nClusters')
    
    % How many neurons can each neuron connect to
    nPostSynCand = nNeuronsPerCluster.*mncVec';
    figure; imagesc(mncVec, nClustersVec, nPostSynCand'); %, 'AlphaData', ~isnan(op'));
    colorbar; set(gca,'YDir','normal')
    title('nPostSynCand'); xlabel('mnc'); ylabel('nClusters')
    
    % What fraction of possible connections are present
    fracPostConn =  (parameters.n_E*parameters.conn_prob) ./nPostSynCand; fracPostConn(fracPostConn>1) = 1;
    figure; imagesc(mncVec, nClustersVec, fracPostConn'); %, 'AlphaData', ~isnan(op'));
    colorbar; set(gca,'YDir','normal')
    title('fracPostConn'); xlabel('mnc'); ylabel('nClusters')
    
    
end




%% Functions

function fcc = fagioloCC(W)
% From Fagiolo, 2007.
% Uses equations 6-8 to calculate the clustering coefficient of a
% directed binary graph W

dtot = (W' + W)* ones( 1, size(W,1) )'; % eq 6
dbi = diag(W^2); % eq 7
nodeCCs = diag([ (W + W')^3]) ./ [2 .* (dtot.*(dtot-1) - 2.*dbi ) ]; % eq 8
fcc = mean(nodeCCs);

end


function [acc, c] = avgClusteringCoefficient(graph)
    %
    % https://www.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/45734/versions/1/previews/cnm/avgClusteringCoefficient.m/index.html
    %
    %Computes the Average Clustering Coefficient for the undirected, unweighted
    %graph input as adjacency matrix 'graph' and returns it in variable 'acc'. 
    %It also returns the local clustering coefficient of each node in 'c'.
    %This implementation does not take into account nodes with zero degree.
    %
    %The definition of clustering coefficient used in this implementation was
    %taken from:
    %
    %   Watts,D.J. and Strogatz,S.H. (1998) Collective dynamics of 
    %   "small-world" networks. Nature, 393, 440-442.
    %
    %INPUT
    %   graph -> The adjacency matrix representation of a graph. It has to be a
    %            NxN matrix where N is the number of nodes in the graph. This
    %            parameter is required.
    %
    %OUTPUT
    %   acc -> Average clustering coefficient of the input graph.
    %   c -> Local clustering coefficient of each of the graph's nodes
    %
    %Example usage:
    %
    %   [acc, c] = avgClusteringCoefficient(my_net);
    %   acc = avgClusteringCoefficient(my_net);
    %   [~, c] = avgClusteringCoefficient(my_net);
    %
    % Copyright (C) Gregorio Alanis-Lobato, 2014

    %% Input parsing and validation
    ip = inputParser;

    %Function handle to make sure the matrix is symmetric
    issymmetric = @(x) all(all(x == x.'));

    %addRequired(ip, 'graph', @(x) isnumeric(x) && issymmetric(x));
    addRequired(ip, 'graph');

    parse(ip, graph);

    %Validated parameter values
    graph = ip.Results.graph;


    %% Clustering coefficient computation

    %Make sure the graph unweighted!!!
    graph(graph ~= 0) = 1; 

    deg = sum(graph, 2); %Determine node degrees
    cn = diag(graph*triu(graph)*graph); %Number of triangles for each node

    %The local clustering coefficient of each node
    c = zeros(size(deg));
    c(deg > 1) = 2 * cn(deg > 1) ./ (deg(deg > 1).*(deg(deg > 1) - 1)); 

    %Average clustering coefficient of the graph
    acc = mean(c(deg > 1));
end

function h = WattsStrogatz(N,K,beta)
% H = WattsStrogatz(N,K,beta) returns a Watts-Strogatz model graph with N
% nodes, N*K edges, mean node degree 2*K, and rewiring probability beta.
%
% beta = 0 is a ring lattice, and beta = 1 is a random graph.

% Connect each node to its K next and previous neighbors. This constructs
% indices for a ring lattice.
s = repelem((1:N)',1,K);
t = s + repmat(1:K,N,1);
t = mod(t-1,N)+1;

% Rewire the target node of each edge with probability beta
for source=1:N    
    switchEdge = rand(K, 1) < beta;
    
    newTargets = rand(N, 1);
    newTargets(source) = 0;
    newTargets(s(t==source)) = 0;
    newTargets(t(source, ~switchEdge)) = 0;
    
    [~, ind] = sort(newTargets, 'descend');
    t(source, switchEdge) = ind(1:nnz(switchEdge));
end

h = graph(s,t);
end

function ci = clusteringcoef(g)
% clusteringcoef    - clustering coefficient of given adjacency matrix
%
%   coefs = clusteringcoef(g) cluster coefficient is ratio of the number of
%   edges Ei between the first neighbors of the vertex i, and the
%   respective number of edges, Ei(max) = ai(ai-1)/2, in the complete graph
%   that can be formed by the nearest neighbors of this vertex:
%
%   g is a graph or an alternatively adjacency matrix.
%
%          2 Ei
%   ci = -----------
%        ai (ai - 1)
%
%   A note that do no have a link with others has clustering coefficient of NaN.
if isa(g, 'graph')
    adj = adjacency(g);
else
    adj = g;
end
n = length(adj);
ci = zeros(1,n);
for k = 1:n
    neighbours = [find(adj(:,k))]';
    neighbours(neighbours==k) = []; % self link deleted
    a = length(neighbours);
    if a == 0; ci(k) = NaN; continue; end
    if a < 2, continue; end
    E = 0;
    for ks = 1:a
        k1 = neighbours(ks);
        for k2 = neighbours(ks+1:a)
            if adj(k1,k2) || adj(k2,k1)
                E = E + 1;
            end
        end
    end
    ci(k) = 2 * E / (a * (a-1));
end
end
