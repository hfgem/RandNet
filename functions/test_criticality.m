function criticality = test_criticality(avalanche_lengths, avalanche_counts)
    %This function is used to test network results for chaotic activity
    %Test 1: Length: # Successive time bins with no loss of activity
    %Test 2: Size: Total # spikes during stretch of time bins with no loss 
    %       of activity
    %Test 3: Average Size / Length: Mean # spikes of a given duration as a 
    %       function of duration
    %The histograms of all of these values can be fit by power laws if the
    %activity is an avalanche: log(xα) = αlog(x)
    %
    %INPUTS:
    %   avalanche_lengths = an array of (1xN) size, where N is the number of
    %                       events of the network spiking, containing the
    %                       length of each avalanche sequence.
    %   avalanche_counts = an array of (1xN) size, where N is the number of
    %                       events of the network spiking, containing the
    %                       number of spikes in each avalanche sequence.
    
    %Grab numbers of events and appropriate histogram bin numbers based on
    %Freedman-Diaconis rule
    N = length(avalanche_lengths); %Number of events
    iqr_len = iqr(avalanche_lengths); %Interquartile range of lengths
    iqr_count = iqr(avalanche_counts); %Interquartile range of counts
    len_bin_wid = 2*iqr_len/(N^(1/3));
    len_bins = [0:len_bin_wid:max(avalanche_lengths)];
    count_bin_wid = 2*iqr_count/(N^(1/3));
    count_bins = [0:count_bin_wid:max(avalanche_counts)];
    
    %Grab average size of avalanche per length stats and get histogram bin
    %numbers as above
    [sort_ava_len, sort_ind] = sort(avalanche_lengths);
    unique_len = unique(sort_ava_len);
    avg_sz_per_len = zeros(size(unique_len));
    for u_l = 1:length(unique_len)
        %Write some code finding the average size per length
        len_ind = sort_ava_len == unique_len(u_l);
        sizes_per_len = avalanche_counts(sort_ind(len_ind));
        avg_sz_per_len(u_l) = mean(sizes_per_len);
    end
    clear u_l len_ind sizes_per_len
    
    %Test for linear relationships
    
    
    %Set the criticality flag based on test
    criticality = 0;
    
    
    

end