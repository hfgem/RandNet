function [D,D_sd] = calculate_trajectory_distances(X)
    %ABOUT: This function takes a structure with multi-dimensional
    %trajectories and calculates, pairwise, their distances in space.
    %
    %INPUTS:
    %   1. X = a structure that is [1,seq] in size, where seq
    %       is the number of sequences. This structure only has 1 field
    %       "sequence", which contains [clusters,t] size matrices where
    %       clusters = the number of clusters (dimensions) and t = the
    %       number of timesteps.
    %OUTPUTS:
    %   1. D = a matrix of average distances calculated between
    %          trajectories, of size [seq,seq], with a diagonal of 0s.
    %   2. D_sd = a matrix of standard deviations of trajectory distance 
    %          vector values, of size [seq,seq], with a diagonal of 0s.
    
    [~,seq] = size(X);
    D = zeros(seq,seq);
    D_sd = zeros(seq,seq);
    for i = 1:seq-1 %Sequence 1
        seq_1 = X(i).sequence;
        [~,t1] = size(seq_1);
        for j = i+1:seq %Sequence 2
            seq_2 = X(j).sequence;
            [~,t2] = size(seq_2);
            t = min([t1,t2]);
            seq_1_short = seq_1(:,1:t);
            seq_2_short = seq_2(:,1:t);
            %Calculate the trajectory distances
            d_vec = sqrt(sum((seq_1_short - seq_2_short).^2,1));
            D(i,j) = mean(d_vec);
            D(j,i) = mean(d_vec);
            D_sd(i,j) = std(d_vec);
            D_sd(j,i) = std(d_vec);
        end    
    end
end