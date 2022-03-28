function [full_ranks, sequence_lengths, nonfiring_neurons, sequence_mat] = create_rank_matrix(spike_sequences, varargin)
%ABOUT: This function takes a spike sequence matrix and generates a matrix
%of ranks from it. Ranks of nonspiking neurons will be ns + 0.5*(n- ns)
%where ns is the number of participating neurons in the sequence and n is
%the total number of neurons in the network.
%
% Treatment of nonspiking neuron ranks can be controlled by optional input
% "nonspike_method"
% middle: ns + 0.5*(n- ns)
% end: ns+1
% endSorted: nonspiking neurons sorted by the index in network
% 0: 0
% nan: nan

    % Defaults for optional parameters
    nonspike_method = 'middle';  % possible options: 'middle', 'end', 'endSorted', 'nan', '0'
    
    % Read in optional parameters, to overwrite above defaults
    for i=1:2:length(varargin)
        switch varargin{i}
            case 'nonspike_method'
                nonspike_method = varargin{i+1};
            otherwise
                error('create_rank_matrix: Unknown input')
        end
    end
    
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
                ns = sum( ~[isnan(seq_r)|seq_r==0] ); % number that are spiking, indicated by a 0 or a nan
                padded_seq = [seq;zeros(n - ns,1)];
                [x,~] = size(seq_r);
                if x == 1
                    seq_r = seq_r';
                end
                
                zero_ind = [isnan(seq_r)|seq_r==0]; % index of neurons that did not spike
                switch nonspike_method
                    case 'middle'
                        new_rank = ns + 0.5*(n - ns);
                    case 'nan'
                        new_rank = nan;
                    case 'end'
                        new_rank = ns+1;
                    case 'endSorted'
                        new_rank = ns+1:n;
                    case '0'
                        new_rank = 0;
                    otherwise 
                        error('Unknown nonspike_method')
                end
                seq_r(zero_ind) = new_rank; %Update ranks of nonfiring neurons to be halfway between seq length and network size
                
                full_ranks = [full_ranks, seq_r]; %#ok<*AGROW>
                sequence_lengths = [sequence_lengths, ns];
                nonfiring_neurons = [nonfiring_neurons, zero_ind];
                sequence_mat = [sequence_mat, padded_seq];
            end % Sequence loop
        end 
    end % Test loop
    
end