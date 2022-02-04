function [short_order_sim, short_order] = order_similarity(sequence_mat, cluster_mat)
    %ABOUT: This function takes a rank matrix and cluster assignment inputs
    %and calculates whether successive neurons are more likely to be from
    %the same cluster as those preceding, or different. This is done in 2
    %ways:
    %   (1) excluding nonspiking neurons summing the amount of cluster
    %   overlap and dividing by the number of pairs for an average overlap
    %   of each sequence
    %   (2) a matrix containing individual pair overlap in short sequences 
    %   (before averaging in (2) above).
    %
    %INPUTS:
    %   1. sequence_mat = a matrix of n x n_s size where n is the total
    %       number of neurons in the network and n_s the number of spike
    %       sequences. Each column contains the order of neurons firing (by
    %       indices)
    %   2. cluster_mat = a matrix of size n_c x n where n_c is the number
    %       of clusters in the network, and n is the number of neurons.
    %       Each row contains a binary marking of whether a neuron (column)
    %       is in a particular cluster (row).
    %
    %OUTPUTS:
    %   1. short_order_sim = vector of size 1 x n_s containing the average
    %       overlap for each short sequence (bullet 2 in about)
    %   2. short_order = matrix of size n x n_s containing the cluster
    %       representation for each pair in each short sequence (bullet 4 in about)

    %Grab relevant parameters
    [n, n_s] = size(sequence_mat);
    
    %Set up variables
    short_order_sim = zeros([1,n_s]);
    short_order = zeros([n,n_s]);
    
    %For each sequence, calculate the cluster similarity in the order of
    %firing
    for i = 1:n_s
        %Pull out the long and short sequences
        full_seq = sequence_mat(:,i);
        short_seq = full_seq(full_seq~= 0);
        %Run through the ordering and calculate the clusters represented in
        %each pair.
        short_total_pairs = length(short_seq) - 1;
        short_overlap = 0;
        short_order_seq = zeros([n,1]);
        for j = 1:short_total_pairs
            c_1 = cluster_mat(:,short_seq(j));
            c_2 = cluster_mat(:,short_seq(j+1));
            pair_overlap = sum(c_1.*c_2); %If the two neurons overlap in more than 1 cluster, that is counted in the sum.
            short_overlap = short_overlap + pair_overlap;
            short_order_seq(j) = pair_overlap;
        end
        short_order_sim(i) = short_overlap/short_total_pairs;
        short_order(:,i) = short_order_seq;
    end    
    
end