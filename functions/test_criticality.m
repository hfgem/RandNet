function criticality = test_criticality(avalanche_lengths, avalanche_counts)
    %This function is used to test network results for chaotic activity
    %Test 1: Length: # Successive time bins with no loss of activity
    %Test 2: Size: Total # spikes during stretch of time bins with no loss 
    %       of activity
    %Test 3: Average Size / Length: Mean # spikes of a given duration as a 
    %       function of duration
    %The histograms of all of these values can be fit by power laws if the
    %activity is an avalanche.
    %We use the MATLAB linear fit function to fit the data and test:
    %   1. Whether the fit is good based on the distribution of the
    %   residuals - if the residuals appear to behave randomly,
    %   it suggests that the model fits the data well. Thus we perform a
    %   KS-Test to determine if the residuals are random-normally
    %   distributed.
    %   2. Whether the slope of the linear fit is sufficiently > 0 and
    %   maintains the relationship to the other slopes needed for
    %   criticality.
    %
    %INPUTS:
    %   avalanche_lengths = an array of (1xN) size, where N is the number of
    %                       events of the network spiking, containing the
    %                       length of each avalanche sequence.
    %   avalanche_counts = an array of (1xN) size, where N is the number of
    %                       events of the network spiking, containing the
    %                       number of spikes in each avalanche sequence.
    
    %Storage for criticality results of 3 tests
    crit_3 = zeros(1,3);
    
    %Grab numbers of events and appropriate histogram bin numbers based on
    %Freedman-Diaconis rule
    N = length(avalanche_lengths); %Number of events
    iqr_len = iqr(avalanche_lengths); %Interquartile range of lengths
    iqr_count = iqr(avalanche_counts); %Interquartile range of counts
    len_bin_wid = 2*iqr_len/(N^(1/3));
    len_bins = (0:len_bin_wid:max(avalanche_lengths));
    len_bin_avgs = (len_bins(1:end-1) + len_bins(2:end))/2;
    [len_hist,~] = histcounts(avalanche_lengths,len_bins);
    count_bin_wid = 2*iqr_count/(N^(1/3));
    count_bins = (0:count_bin_wid:max(avalanche_counts));
    count_bin_avgs = (count_bins(1:end-1) + count_bins(2:end))/2;
    [count_hist,~] = histcounts(avalanche_counts,count_bins);
    
    %First Test: Power law of #avalanches of given duration as a function
    %           of duration
    len_dataset = table(log(len_bin_avgs'),log(len_hist'));
    len_mdl = fitlm(len_dataset);
    resid = len_mdl.Residuals.Standardized;
    h = kstest(resid); %Test 1 that residuals are standard-normally distributed
    len_power = 1;
    if h ~= 0
        len_power = 0;
    end
    len_slope = -1*len_mdl.Coefficients.Estimate(2);
    if abs(len_slope) < 0.5 %Test 2 that slope is large enough
        len_power = 0;
    end
    crit_3(1) = len_power;
    
    %Second Test: Power law of #avalanches of a given size as a function of
    %           size
    count_dataset = table(log(count_bin_avgs'),log(count_hist'));
    count_mdl = fitlm(count_dataset);
    resid = count_mdl.Residuals.Standardized;
    h = kstest(resid); %Test 1 that residuals are standard-normally distributed
    count_power = 1;
    if h ~= 0
        count_power = 0;
    end
    count_slope = -1*count_mdl.Coefficients.Estimate(2);
    if abs(count_slope) < 0.5 %Test 2 that slope is large enough
        count_power = 0;
    end
    crit_3(2) = count_power;
    
    %Third Test: Grab average size of avalanche per length stats and fit a 
    %           line to it
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
    avg_sz_dataset = table(log(unique_len'),log(avg_sz_per_len'));
    avg_sz_mdl = fitlm(avg_sz_dataset);
    resid = avg_sz_mdl.Residuals.Standardized;
    h = kstest(resid); %Test 1 that residuals are standard-normally distributed
    avg_sz_power = 1;
    if h ~= 0
        %The residuals are not normally distributed and so the linear fit
        %is good enough
        avg_sz_power = 0;
    end
    %       Test 2 that, if avalanches are scale-free, then the exponents are not 
    %       independent, but satisfy the equation (slope_1 - 1)/(slope_2 -
    %       1) = 1/slope_3
    avg_sz_slope = avg_sz_mdl.Coefficients.Estimate(2);
    slope_ratio = (len_slope - 1)/(count_slope - 1);
    slope_dist = slope_ratio - 1/avg_sz_slope;
    if slope_dist > 0.05*slope_ratio
        avg_sz_power = 0;
    end
    
    %Set the criticality flag based on three tests
    criticality = prod(crit_3,'all');

end