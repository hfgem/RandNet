function [network] = create_clusters(parameters, varargin)
    %_________
    %ABOUT: This function generates the network clusters and connections
    %based on the number of neurons, number of clusters, number of neurons
    %per cluster, and the probability of connections within a cluster.
    %
    %INPUTS:
    %   parameters = struct containing the following
    %       n = number of neurons in the network
    %       clusters = number of clusters of neurons
    %       cluster_n = number of neurons in a cluster
    %       cluster_prob = intra-cluster connection probability
    %   seed = random number generator seed
    %   include_all = binary value of whether to include all neurons in
    %       clusters (1) or not (0).
    %   global_inhib = binary value of whether to include a global
    %       inhibition to the connectivity matrix.
    %OUTPUTS:
    %   cluster_mat = A binary [clusters x n] matrix of which neurons are in
    %                 which cluster
    %   conns = An [n x n] matrix of which neurons are connected to each
    %                 other, with values greater than 1 implying stronger
    %                 connectivity
    %_________

        
    %% Defaults for optional parameters
    seed =  randi(10^6); % rng seed, to replicate network construction
    include_all = 1;    % if 1, all neurons are included in at least 1 cluster by reducing high participating cells
    global_inhib = 1;   % if 1, I-cells have widespread, unclustered connectivity
   
    
    %% Read in optional parameters, to overwrite above defaults
    for i=1:2:length(varargin)
        switch varargin{i}
            case 'seed'
                seed = varargin{i+1};
            case 'include_all'
                include_all = varargin{i+1};
            case  'global_inhib'
                global_inhib = varargin{i+1};
            otherwise
                error('create_clusters: Unknown input')
        end
    end

    
    %% Main:
    rng(seed)
    
    %Decide which neurons are inhib and excit 
    all_indices = [1:parameters.n];
    I_indices = datasample(all_indices,parameters.n_I,'Replace',false); %indices of inhibitory neurons
    E_indices = find(~ismember(all_indices,I_indices)); %indices of excitatory neurons

    %Create clusters by randomly selecting cluster_n neurons for each
    cluster_mat = zeros(parameters.clusters,parameters.n);
    for i = 1:parameters.clusters %set clusters
        ord = randperm(parameters.n,parameters.cluster_n); %random connectivity
        cluster_mat(i,ord) = 1; %note which neurons are in a cluster
    end
    clear ord i
    
    if include_all == 1 
        %Add back in those neurons that are not placed in a cluster, by
        %removing a place from another neuron with a high presence - this
        %section can be removed if you'd like some neurons to be unconnected
        ind_non = find(sum(cluster_mat) == 0);
        ind_high = find(sum(cluster_mat) > 2);
        for i = ind_non
            clust_place = randi(parameters.clusters);
            ind_inclust = find(cluster_mat(clust_place,:));
            try %#ok<TRYNC>
                val_steal = datasample(intersect(ind_inclust,ind_high),1);
                cluster_mat(clust_place,i) = 1;
                cluster_mat(clust_place,val_steal) = 0;
                ind_high = find(sum(cluster_mat) > 2);
            end
        end
    elseif include_all==2
        % Same as include_all==1, but ind_high threshold is changed to 1
        ind_non = find(sum(cluster_mat) == 0);
        ind_high = find(sum(cluster_mat) > 1);
        for i = ind_non
            clust_place = randi(parameters.clusters);
            ind_inclust = find(cluster_mat(clust_place,:));
            try %#ok<TRYNC>
                val_steal = datasample(intersect(ind_inclust,ind_high),1);
                cluster_mat(clust_place,i) = 1;
                cluster_mat(clust_place,val_steal) = 0;
                ind_high = find(sum(cluster_mat) > 1);
            end
        end
    end
    
    %Find the matrix of total connectivity - it will have integer values of
    %1 or more, representing the strength of connection between neurons
    conns = zeros(parameters.n);
    for i = 1:parameters.clusters
        ord = find(cluster_mat(i,:));
        ord_len = length(ord);
        conns(ord,ord) = conns(ord,ord) + (rand(ord_len,ord_len) <= parameters.cluster_prob);
    end
    
    %Remove self-connectivity
    for i = 1:parameters.n
        conns(i,i) = 0;
    end
    
    % If global_inhib, overwrite all inhibitory outputs+inputs
    if global_inhib
        conns(I_indices,:) = (rand([parameters.n*(1-parameters.p_E), parameters.n]) < parameters.p_I);
        conns(:,I_indices) = (rand(parameters.n, [parameters.n*(1-parameters.p_E)]) < parameters.p_I);
    end
    
    % Input strengths
    input1 = lognrnd(parameters.W_gin, parameters.cueSigma, parameters.n, 1); % location cue 1 strengths
	input2 = lognrnd(parameters.W_gin, parameters.cueSigma, parameters.n, 1); % location cue 2 strengths
    contextInput = lognrnd(parameters.W_gin, parameters.cueSigma/4, parameters.n, 1); % context cue strength
    % contextInput = parameters.Win_mean+(sqrt(parameters.Win_var)*randn(parameters.n, 1)); % context cue strength
    %contextInput = lognrnd(log(parameters.W_gin), parameters.cueSigma, parameters.n, pfsim.nEnvironments);
    
    % I cells don't receive spatially modulated input
    input1(I_indices) = 0;
    input2(I_indices) = 0;
    
    %SAVE NETWORK STRUCTURE
    network = struct;
    network.cluster_mat = cluster_mat;
    network.conns = conns;
    network.seed = seed;
    network.n_EE = sum(network.conns(E_indices,E_indices),'all'); %number of E-E connections
    network.n_EI = sum(network.conns(E_indices,I_indices),'all'); %number of E-I connections
    network.n_II = sum(network.conns(I_indices,I_indices),'all'); %number of I-I connections
    network.n_IE = sum(network.conns(I_indices,E_indices),'all'); %number of I-E connections
    network.all_indices = all_indices;
    network.I_indices = I_indices;
    network.E_indices = E_indices;
    
    network.spatialInput{1} = input1;
    network.spatialInput{2} = input2;
    network.contextInput = contextInput;
    
end