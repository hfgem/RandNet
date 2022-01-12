function [full_ranks, sequence_lengths, nonfiring_neurons] = create_rank_matrix(spike_sequences)
%ABOUT: This function takes a spike sequence matrix and generates a matrix
%of ranks from it. Ranks of nonspiking neurons will be ns + 0.5*(n- ns)
%where ns is the number of participating neurons in the sequence and n is
%the total number of neurons in the network.
    full_ranks = [];
    sequence_lengths = [];
    nonfiring_neurons = [];
    for i = 1:length(spike_sequences)
        if ~isempty(spike_sequences(i).spike_order)
            fieldnames = fields(spike_sequences(i).spike_ranks);
            for j = 1:length(fieldnames)
                seq = spike_sequences(i).spike_ranks.(fieldnames{j});
                n = length(seq);
                ns = sum(seq~=0);
                [x,~] = size(seq);
                if x == 1
                    seq = seq';
                end
                new_rank = ns + 0.5*(n - ns);
                zero_ind = seq == 0;
                seq(zero_ind) = new_rank; %Update ranks of nonfiring neurons to be halfway between seq length and network size
                full_ranks = [full_ranks, seq]; %#ok<*AGROW>
                sequence_lengths = [sequence_lengths, ns];
                nonfiring_neurons = [nonfiring_neurons, zero_ind];
            end    
        end
    end  
end