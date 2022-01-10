function S = generate_shuffled_cluster_trajectories(X,n)
    %ABOUT: This function takes a structure with multi-dimensional
    %trajectories and generates n shuffled trajectories
    %
    %INPUTS:
    %   1. X = a structure that is [1,seq] in size, where seq
    %       is the number of sequences. This structure only has 1 field
    %       "sequence", which contains [clusters,t] size matrices where
    %       clusters = the number of clusters (dimensions) and t = the
    %       number of timesteps.
    %   2. n = the number of shuffles to generate
    %OUTPUTS:
    %   1. S = a structure that is [1,n] in size. This structure only has 1
    %   field "sequence", which contains [clusters,t] size matrices where
    %   clusters = the number of clusters (dimensions) and t = the number
    %   of timesteps.
    
    [~,seq] = size(X);
    S = struct;
    for i = 1:n
        samp_i = randsample(seq,1,false); %Randomly sample a sequence to shuffle
        true_seq = X(samp_i).sequence; %Grab true sequence
        [clusters,seq_length] = size(true_seq); %Grab dimensions of sequence
        new_seq = zeros(size(true_seq)); %Storage of shuffled sequence
        for j = 1:seq_length %For each timestep, shuffle the true sequence
            shuffle = randsample(clusters,clusters,false);
            new_seq(:,j) = true_seq(shuffle,j);
        end    
        S(i).sequence = new_seq;
        clear samp_i true_seq clusters seq_length new_seq j shuffle
    end    
end