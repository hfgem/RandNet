function [cluster_mat, conns] = create_clusters(parameters, seed, include_all)
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
    %OUTPUTS:
    %   cluster_mat = A binary [clusters x n] matrix of which neurons are in
    %                 which cluster
    %   conns = An [n x n] matrix of which neurons are connected to each
    %                 other, with values greater than 1 implying stronger
    %                 connectivity
    %_________

    rng(seed)
    
    %Create clusters by randomly selecting cluster_n neurons for each
    cluster_mat = zeros(parameters.clusters,parameters.n);
    for i = 1:parameters.clusters %set clusters
        ord = randperm(parameters.n,parameters.cluster_n); %random connectivity
        cluster_mat(i,ord) = 1; %note which neurons are in a cluster
    end
    clear ord i

    if include_all
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
end