function test_criticality_hpcc(avalanche_lengths, avalanche_counts, parameters, ithParamSet)
    %This function is used to plot criticality results for network activity
    %Test 1: Length: # Successive time bins with no loss of activity
    %Test 2: Size: Total # spikes during stretch of time bins with no loss 
    %       of activity
    %Test 3: Average Size / Length: Mean # spikes of a given duration as a 
    %       function of duration
    %The histograms of all of these values can be fit by power laws if the
    %activity is an avalanche. This function simply plots the results.
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
    if N > 10
        iqr_len = iqr(avalanche_lengths); %Interquartile range of lengths
        iqr_count = iqr(avalanche_counts); %Interquartile range of counts
        len_bin_wid = 2*iqr_len/(N^(1/3));
        len_bins = (0:len_bin_wid:max(avalanche_lengths));
        if length(len_bins) < 5
           len_bins = linspace(0,max(avalanche_lengths),10);
        end    
        len_bin_avgs = (len_bins(1:end-1) + len_bins(2:end))/2;
        [len_hist,~] = histcounts(avalanche_lengths,len_bins);
        len_nonzero = find(len_hist);
        count_bin_wid = 2*iqr_count/(N^(1/3));
        count_bins = (0:count_bin_wid:max(avalanche_counts));
        if length(count_bins) < 5
           count_bins = linspace(0,max(avalanche_lengths),10);
        end    
        count_bin_avgs = (count_bins(1:end-1) + count_bins(2:end))/2;
        [count_hist,~] = histcounts(avalanche_counts,count_bins);
        count_nonzero = find(count_hist);

        if length(len_nonzero) && length(count_nonzero) > 1 %#ok<ISMT>
            [sort_ava_len, sort_ind] = sort(avalanche_lengths);
            unique_len = unique(sort_ava_len);
            avg_sz_per_len = zeros(size(unique_len));
            for u_l = 1:length(unique_len)
                %Code finding the average size per length
                len_ind = sort_ava_len == unique_len(u_l);
                sizes_per_len = avalanche_counts(sort_ind(len_ind));
                avg_sz_per_len(u_l) = mean(sizes_per_len);
            end
            nonzero_avg_sz = find(avg_sz_per_len.*unique_len);
            clear u_l len_ind sizes_per_len

            %Plot Results
            if parameters.plotResults == 1
                f = figure;
                subplot(1,3,1) %Test 1
                scatter(log(len_bin_avgs(len_nonzero)),log(len_hist(len_nonzero)))
                xlabel('Log(Length)')
                ylabel('Log(Length Counts)')

                subplot(1,3,2)
                scatter(log(count_bin_avgs(count_nonzero)),log(count_hist(count_nonzero)))
                xlabel('Log(Size)')
                ylabel('Log(Size Counts)')

                subplot(1,3,3)
                scatter(log(unique_len(nonzero_avg_sz)),log(avg_sz_per_len(nonzero_avg_sz)))
                xlabel('Log(Length)')
                ylabel('Log(Average Size Per Length)')

                f.Units = 'Normalized';
                f.OuterPosition = [0 0 1 1];
                if ~isfolder(strcat(parameters.save_path,'/crit_plots'))
                    mkdir(parameters.save_path,'/crit_plots')
                end
                savefig(f,strcat(parameters.save_path,'/crit_plots','/sim_',string(ithParamSet),'_plaw.fig'))
                saveas(f,strcat(parameters.save_path,'/crit_plots','/sim_',string(ithParamSet),'_plaw.jpg'))
                close(f)
            end 
        end
    end
end