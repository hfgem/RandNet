function [matching_index, matching_index_mod] = calculate_trajectory_similarity_mi(n, ...
    viable_inits, network_spike_sequences, penalize_nonspike)
    %_________
    %ABOUT: This function calculates the Matching Index (MI) from Vas et al
    %for pairs of sequences including and excluding nonfiring neurons.
    %
    %INPUTS:
    %   n = number of neurons in the simulation
    %   viable_inits = vector of indices of viable initializations in the 
    %           simulation data
    %   network_spike_sequences = a [1 x inits] size structure including 
    %           all of the initializations, their events, the spike order
    %           for each event, the spike ranks for each event, and the
    %           nonspiking neurons
    %   penalize_nonspike = a binary flag. If = 1, it will act on (1) the
    %           matching index calculation including nonspiking neurons by
    %           making sequence pairs where a neuron is nonspiking in one
    %           sequence, and spiking in another, dissimilar, and (2) the
    %           matching index mod calculation excluding nonspiking neurons
    %           by making sequence pairs where a neuron is nonspiking in
    %           one sequence and spiking in another, dissimilar.
    %
    %OUTPUTS:
    %   matching_index = a [num_seq x num_seq] matrix of sequence pair
    %           matching indices calculated for sequences including
    %           nonspiking neurons
    %   matching_index_mod = a [num_seq x num_seq] matrix of sequence pair
    %           matching indices calculated for sequences excluding
    %           nonspiking neurons
    %_________
    
    %First store each event's order information into a new matrix
    all_orders = [];
    all_nonspiking = [];
    for i = viable_inits
        if ~isempty(network_spike_sequences(i).spike_order) %second check
            sequences = network_spike_sequences(i).spike_order(1);
            sequence_names = fieldnames(sequences);
            [num_seq,~] = size(sequence_names);
            for j = 1:num_seq
                order_vals = network_spike_sequences(i).spike_order.(sequence_names{j});
                %now find the nonspiking neurons and add them to the end
                nonspiking_neurons = find(network_spike_sequences(i).nonspiking_neurons.(sequence_names{j}));
                new_order_vals = [order_vals', nonspiking_neurons];
                all_orders = [all_orders; new_order_vals]; %#ok<AGROW>
                all_nonspiking = [all_nonspiking; network_spike_sequences(i).nonspiking_neurons.(sequence_names{j})'];  %#ok<AGROW>
            end  
        end
    end
    clear i sequences sequence_names num_seq j rank_vals
    [num_seq, ~] = size(all_orders);
    all_orders = all_orders'; %rotate for corr function
    
    %======MATCHING INDICES - INCLUDING NONSPIKING======
    
    matching_index = zeros(num_seq,num_seq);
    for s1 = 1:num_seq %for each pair of sequences calculate MI
        seq1 = all_orders(:,s1);
        for s2 = 1:num_seq
            seq2 = all_orders(:,s2);
            match = 0;
            non_match = 0;
            for i = 1:n %first neuron index
                s1_i = find(seq1 == i);
                s2_i = find(seq2 == i);
                s1_i_n = all_nonspiking(i,s1);
                s2_i_n = all_nonspiking(i,s2);
                for j = 1:n %second neuron index
                    s1_j = find(seq1 == j);
                    s2_j = find(seq2 == j);
                    s1_j_n = all_nonspiking(j,s1);
                    s2_j_n = all_nonspiking(j,s2);
                    %now compare order
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
                end
            end
            MI = (match - non_match)/(match + non_match);
            matching_index(s1,s2) = MI;
            matching_index(s2,s1) = MI;
        end
    end
    %clear s1 seq1 s2 seq2 match non_match i s1_i s2_i j s1_j s2_j o1 o2 MI
    
    %======MATCHING INDICES - NOT INCLUDING NONSPIKING======

    matching_index_mod = zeros(num_seq,num_seq);
    for s1 = 1:length(num_seq) %for each pair of sequences calculate MI
        seq1 = all_orders(:,s1);
        seq1_nonspike = find(all_nonspiking(:,s1) == 1);
        for s2 = 1:length(num_seq)
            seq2 = all_orders(:,s2);
            seq2_nonspike = find(all_nonspiking(:,s2) == 1);
            %find all nonspiking neurons between the two sets and remove
            %from both spike orders for comparison
            if penalize_nonspike == 1
                seq1 = setdiff(seq1, seq2_nonspike','stable'); %Keeps nonspiking neurons unique to seq1 in the sequence
                seq2 = setdiff(seq2, seq1_nonspike','stable'); %Keeps nonspiking neurons unique to seq2 in the sequence
            else
                nonspiking_total = unique([seq2_nonspike',seq1_nonspike']);
                seq1 = setdiff(seq1, nonspiking_total','stable');
                seq2 = setdiff(seq2, nonspiking_total','stable');
            end
            match = 0;
            non_match = 0;
            %if penalize_nonspike == 1, if a neuron is in one sequence but
            %not another, the find function will return [], and the
            %comparisons will give a value of 1x0 logical. This 1x0 logical
            %compared with anything, including a 1x0 logical, will result
            %in non_match + 1;
            for i = 1:n %first neuron index
                s1_i = find(seq1 == i);
                s2_i = find(seq2 == i);
                for j = 1:n %second neuron index
                    if j ~= i
                        s1_j = find(seq1 == j);
                        s2_j = find(seq2 == j);
                        %now compare order
                        o1 = s1_i > s1_j;
                        o2 = s2_i > s2_j;
                        if o1 == o2
                            match = match + 1;
                        else
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
    clear s1 seq1 seq1_nonspike s2 seq2 seq2_nonspike nonspiking match ...
        non_match i s1_i s2_i j s1_j s2_j o1 o2 MI
    
end