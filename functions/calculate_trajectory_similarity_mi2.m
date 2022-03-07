function [matching_index, matching_index_mod] = calculate_trajectory_similarity_mi2(X, ...
    X_n, penalize_nonspike)
    %________
    %ABOUT: This function takes in a matrix of sequences and calculates the
    %   Matching Index (MI) from Vas et al for pairs of sequences including
    %   and excluding nonfiring neurons.
    %
    %INPUTS:
    %   1. X = a matrix that is n x m in size, where n is the number of
    %           neurons, and m is the number of sequences to compare. Each
    %           (i,j) value is the rank of the ith neuron in sequence j. If
    %           a neuron did not spike in a sequence, the rank will be
    %           halfway between the total number of spiking neurons (n_s)
    %           and the total number n (n_s + (n-n_s)/2).
    %   2. X_n = a matrix with the same dimensions as X, but containing
    %           binary values of whether (0) or not (1) a neuron spiked.
    %   3. penalize_nonspike = a binary flag. If = 1, it will penalize
    %       nonspiking neurons 
    
    [n, num_seq] = size(X);
    
    %======MATCHING INDICES - INCLUDING NONSPIKING======
    matching_index = zeros(num_seq, num_seq);
    for s1 = 1:num_seq %for each pair of sequences calculate MI
        seq1 = X(:,s1);
        for s2 = 1:num_seq
            seq2 = X(:,s2);
            match = 0;
            non_match = 0;
            for i = 1:n %first neuron index
                s1_i = seq1(i);
                s2_i = seq2(i);
                s1_i_n = X_n(i,s1);
                s2_i_n = X_n(i,s2);
                for j = 1:n %second neuron index
                    s1_j = seq1(j);
                    s2_j = seq2(j);
                    s1_j_n = X_n(j,s1);
                    s2_j_n = X_n(j,s2);
                    %now compare order
                    try
                        o1 = s1_i > s1_j;
                        o2 = s2_i > s2_j;
                        if penalize_nonspike == 1
                            if s1_i_n == s2_i_n && s1_j_n == s2_j_n
                                if o1 == o2
                                    match = match + 1;
                                else
                                    non_match = non_match + 1;
                                end
                            else
                                non_match = non_match + 1;
                            end
                        else
                            if o1 == o2
                                match = match + 1;
                            else
                                non_match = non_match + 1;
                            end
                        end
                    catch
                        non_match = non_match + 1;
                    end    
                end
            end
            MI = (match - non_match)/(match + non_match);
            matching_index(s1,s2) = MI;
            matching_index(s2,s1) = MI;
        end
    end
    
    
    %======MATCHING INDICES - NOT INCLUDING NONSPIKING======
    matching_index_mod = zeros(num_seq, num_seq);
    for s1 = 1:num_seq %for each pair of sequences calculate MI
        seq1 = X(:,s1);
        seq1_nonspike = X_n(:,s1);
        for s2 = 1:num_seq
            seq2 = X(:,s2);
            seq2_nonspike = X_n(:,s2);
            %binary vector forms of nonspiking neurons
            overlap_nonspike = seq1_nonspike.*seq2_nonspike;
            nonspiking_total = (seq1_nonspike + seq2_nonspike) > 0;
            %storage for match / non-match counts
            match = 0;
            non_match = 0;
            %if penalize_nonspike == 1, if a neuron is in one sequence but
            %not another, the find function will return [], and the
            %comparisons will give a value of 1x0 logical. This 1x0 logical
            %compared with anything, including a 1x0 logical, will result
            %in non_match + 1;
            for i = 1:n %first neuron index ranks
                s1_i = seq1(i);
                s2_i = seq2(i);
                for j = 1:n %second neuron index ranks
                    if j ~= i
                        s1_j = seq1(j);
                        s2_j = seq2(j);
                        %now compare order
                        try
                            o1 = s1_i > s1_j;
                            o2 = s2_i > s2_j;
                            if penalize_nonspike == 1 %Excludes overlapping nonspiking neurons across both sequences
                                i_n = overlap_nonspike(i);
                                j_n = overlap_nonspike(j);
                                if i_n == 0 && j_n == 0 %both are not in the overlapping population
                                    if o1 == o2
                                        match = match + 1;
                                    else
                                        non_match = non_match + 1;
                                    end
                                end
                            else %Excludes all nonspiking neurons across both sequences
                                i_n = nonspiking_total(i);
                                j_n = nonspiking_total(j);
                                if i_n == 0 && j_n == 0 %both are not nonspiking
                                    if o1 == o2
                                        match = match + 1;
                                    else
                                        non_match = non_match + 1;
                                    end
                                end
                            end
                        catch
                            non_match = non_match + 1;
                        end
                    end
                end
            end
            MI = (match - non_match)/(match + non_match);
            matching_index_mod(s1,s2) = MI;
            matching_index_mod(s2,s1) = MI;
        end
    end
    
end