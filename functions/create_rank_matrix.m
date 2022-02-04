function [full_ranks, sequence_lengths, nonfiring_neurons, sequence_mat] = create_rank_matrix(spike_sequences)
%ABOUT: This function takes a spike sequence matrix and generates a matrix
%of ranks from it. Ranks of nonspiking neurons will be ns + 0.5*(n- ns)
%where ns is the number of participating neurons in the sequence and n is
%the total number of neurons in the network.
    sequence_mat = [];
    full_ranks = [];
    sequence_lengths = [];
    nonfiring_neurons = [];
    for i = 1:length(spike_sequences)
        if ~isempty(spike_sequences(i).spike_order)
            fieldnames = fields(spike_sequences(i).spike_ranks);
            for j = 1:length(fieldnames)
                seq_r = spike_sequences(i).spike_ranks.(fieldnames{j});
                seq = spike_sequences(i).spike_order.(fieldnames{j});
                n = length(seq_r);
                ns = sum(seq_r~=0);
                padded_seq = [seq;zeros(n - ns,1)];
                [x,~] = size(seq_r);
                if x == 1
                    seq_r = seq_r';
                end
                new_rank = ns + 0.5*(n - ns);
                zero_ind = seq_r == 0;
                seq_r(zero_ind) = new_rank; %Update ranks of nonfiring neurons to be halfway between seq length and network size
                full_ranks = [full_ranks, seq_r]; %#ok<*AGROW>
                sequence_lengths = [sequence_lengths, ns];
                nonfiring_neurons = [nonfiring_neurons, zero_ind];
                sequence_mat = [sequence_mat, padded_seq];
            end    
        end
    end  
end