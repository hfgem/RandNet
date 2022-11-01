
%% Plot combined place fields at a particular parameter point
% load('C:\Users\Jordan\Box\Data\Replay project\RandNet\temp, new PF sim code grids\results_2022-07-12T06-58.mat')
% load('C:\Users\Jordan\Box\Data\Replay project\RandNet\temp, new PF sim code grids\results_2022-07-21T15-05.mat')
% load('../results/randnet_PF_param_tests/results_2022-07-21T15-05.mat')

if ispc
    addpath('C:\Users\Jordan\Box\Data\RandNet-Data\temp, new PF sim code grids')
elseif ismac
    addpath('/Users/jordan/Library/CloudStorage/Box-Box/Data/RandNet-Data/temp, new PF sim code grids')
else
    disp('error')
end


% load('results_2022-07-22T18-29.mat') % Primary clustersXmnc grid (smaller mnc values)
% load('results_2022-07-13T11-44.mat') % clustersXmnc grid up to (5,5)

if isfolder('functions')
    addpath('functions')
else
    addpath( fullfile( '..', 'functions' ) )
end

% paramSetInds = [3, 3; 3 4; 4,2]
% paramSetInds = combvec([1:size(resultsStruct, 1)], [1:size(resultsStruct, 2)])'
variedParam(:).name
variedParam(:).range
paramSetInds = combvec([2], [4])' % best example, for 2022-07-22T18-29 data
% paramSetInds = combvec([1], [1])' % stationary example, for 2022-07-22T18-29 data
% paramSetInds = combvec([9], [1])' % diffuse example, for 2022-07-22T18-29 data

parameters.n
parameters.del_G_sra
variedParam(1).name
variedParam(2).name

% Analysis and plotting parameters
minPeakRate = 2; % minimum peak PF rate to be considered a place cell
useWeightedDecode = 0 % slope and R2 for correlations by either peak prob or weighted prob
nEventsToPlot = 12; 
useMeanPFDensity = true

% Analysis and plotting options
calcScore = false
plotExtraPlots = false % if 1, plot place fields of every network
plotNetsBest = true
plotPvalMat = false 
plotClusterPFs = false
plotNetStruct = false
plotPFStats = false

%% Loop over parameter sets
rng(1); tic
for ithParamSet = 1:size(paramSetInds, 1)
    
    ithParamSet
	ithParam1 = paramSetInds(ithParamSet,1);
    disp(['Param1=', num2str(ithParam1), ': ', variedParam(1).name, '=', num2str(variedParam(1).range(ithParam1))])
	ithParam2 = paramSetInds(ithParamSet,2);
    disp(['Param2=', num2str(ithParam2), ': ', variedParam(2).name, '=', num2str(variedParam(2).range(ithParam2))])
        
    for i = 1:size(variedParam, 2)
    	parameters.(variedParam(i).name) = variedParam(i).range(paramSetInds(i));
    end
        
    temp = nan(1, num_nets);
    
    PFmatE_all = [];
    netInd_all = []; % Index of which net each cell came from
    nClustMemb_all = []; % number of clusters each cell is a member of
    PFscores_all = [];
    bestEventpVals = inf(nEventsToPlot, 1); % nEventsToPlot best event pvals
    bestEventAbsrVals = inf(nEventsToPlot, 1); % nEventsToPlot best event pvals
    bestEventjdVals = inf(nEventsToPlot, 1); % nEventsToPlot best event pvals
    bestEvents_decode = cell(nEventsToPlot, 1);
    bestEvents_relRank = cell(nEventsToPlot, 1);
    bestEvents_relRankColor = cell(nEventsToPlot, 1);
    nEvents = 0;
    cluster_matAll = [];
    allEventRs = [];
    allEventMaxJumps = [];
    allShuffleRs = [];
    allShuffleMaxJumps = [];
    allEventLengths = [];
    allShuffleLengths = [];
    avgDecodePos = zeros(size(pfsim.gridxvals));
    for ithNet = 1:size(resultsStruct, 3)

        % Get matrix of PFs
        PFmat = PFresultsStruct(ithParam1, ithParam2, ithNet).results{1}.linfields;
        E_indices = PFresultsStruct(ithParam1, ithParam2, ithNet).results{1}.E_indices;
        netSeed = PFresultsStruct(ithParam1, ithParam2, ithNet).results{1}.netSeed;

        PFmatE_all = [PFmatE_all; PFmat(E_indices,:)];
        netInd_all = [netInd_all; repmat(ithNet, numel(E_indices), 1)];

        % Recreate network
        netParams=parameters;
        for i = 1:size(variedParam, 2)
            netParams.(variedParam(i).name) = variedParam(i).range(paramSetInds(i));
        end
        netParams = set_depedent_parameters(netParams);
        network = create_clusters(netParams, 'seed', netSeed, 'include_all', netParams.include_all, 'global_inhib', netParams.global_inhib);
        if ~all(network.E_indices==E_indices); error('Incorrect network'); end
        
        nClustMemb_all = [nClustMemb_all; sum(network.cluster_mat(:,E_indices), 1)'];
        cluster_matAll = [cluster_matAll; network.cluster_mat(:,E_indices)'] ;

        tBinSz = resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.tBinSz;
        
        
        %% Calculate all pair-wise place field distances and their connectivity
        % !!! TODO !!! 
        % Note: need to normalize by density of Place-field-peak-distances
        % TODO: add option for peak PF vs expected PF pos
        % TODO: validate that connection probabilities are correct with parameterization
        myPlotSettings(4, 3, 2, 14, [], [], 2) % ppt format
        
        nValidPairs = 0;
        posDepConProb = zeros(size(pfsim.gridxvals)); % Position depenedent connection probability
        rand_posDepConProb = zeros(size(pfsim.gridxvals)); % Rand-pos Position depenedent connection probability
        rand2_posDepConProb = zeros(size(pfsim.gridxvals)); % Rand-pos Position depenedent connection probability
        for ithCell = 1:numel(E_indices)
            for jthCell = (ithCell+1):numel(E_indices)
                PFi = PFmatE_all(ithCell,:);
                PFj =  PFmatE_all(jthCell,:);
                
                if max(PFi)<minPeakRate || max(PFj)<minPeakRate
                    continue
                end
                nValidPairs = nValidPairs+1;
                [iMax, iInd] = max(PFi);
                [jMax, jInd] = max(PFj);
                pairDist = abs(iInd-jInd)+1;
                pairConns = network.conns(E_indices(ithCell), E_indices(jthCell)) + ...
                             network.conns(E_indices(jthCell), E_indices(ithCell));
                posDepConProb(pairDist) = posDepConProb(pairDist)+ pairConns;
                rand_posDepConProb(randi(numel(rand_posDepConProb))) = rand_posDepConProb(pairDist)+ pairConns;
                rand2_posDepConProb(pairDist) = rand2_posDepConProb(pairDist) + [rand(1)>0.5];
            end
        end
        posDepConProb = posDepConProb./nValidPairs;
        rand_posDepConProb = rand_posDepConProb./nValidPairs;
        rand2_posDepConProb = rand2_posDepConProb .* [sum(posDepConProb)/sum(rand2_posDepConProb)];
        
        xPosVals = [0:(numel(rand_posDepConProb)-1)];
        figure; hold on; title(['Network #', num2str(ithNet)])
        plot(xPosVals, rand_posDepConProb, 'k:'); plot(xPosVals, rand2_posDepConProb, 'b--'); 
        plot(xPosVals, posDepConProb, 'r');
        legend({'Shuffle peak-dist.', 'Shuffle conn.-prob.', 'Actual'}, 'Location', 'Best')
        xlabel('PF peak distance (2 cm bins)'); ylabel('Connection prob.');
        drawnow
        
        %keyboard
        
        
        %%

        %Accumulate best events, if there are any
        if ~isempty(resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.pvalue)
            

            if useWeightedDecode
                rsqrs_shuffle = vertcat(resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.shuffle_weightedR2{:}); 
                rsqrs_shuffle = rsqrs_shuffle(:,1);
                rsqrs_shuffle = cellfun(@transpose,rsqrs_shuffle,'UniformOutput',false); 
                rsqrs_preplay = resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.weightedR2(:,1); 
                
                % pvals_preplay = resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.pvalue(:,1); % pvalue, relative to shuffle
                pvals_preplay = nan(size(rsqrs_preplay));
                for ithEvent=1:numel(rsqrs_preplay)
                    pvals_preplay(ithEvent) = mean(rsqrs_preplay(ithEvent)<rsqrs_shuffle{ithEvent});
                end
            else
                rsqrs_shuffle = vertcat(resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.shuffle_rsquare{:}); 
                rsqrs_shuffle = rsqrs_shuffle(:,1);
                rsqrs_preplay = resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.rsquare(:,1); 
                pvals_preplay = resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.pvalue(:,1); % pvalue, relative to shuffle
            end
            
            maxJumps_preplay = resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.maxJump(:,1); 
            maxJumps_shuffle = vertcat(resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.shuffle_maxJump{:}); 
            maxJumps_shuffle = maxJumps_shuffle(:,1);
            allEventRs = [allEventRs; rsqrs_preplay];
            allEventMaxJumps = [allEventMaxJumps; maxJumps_preplay];
            allShuffleRs = [allShuffleRs; rsqrs_shuffle];
            allShuffleMaxJumps = [allShuffleMaxJumps; maxJumps_shuffle];
            
            % Extract nTimeBins for each event
            allEventpMats = vertcat(resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.pMat{:});
            firstTrajpMats = vertcat(allEventpMats{:,1});
            firstTrajpMatLengths = arrayfun(@(x) numel(x.timevec), firstTrajpMats);
            allEventLengths = [allEventLengths; firstTrajpMatLengths];
            
            allShuffleLengths = [allShuffleLengths; repelem(firstTrajpMatLengths, numel(rsqrs_shuffle{1}))]; % Note: Shuffles are same length as the corresponding original event
            
            % Accumulate decode probability across spatial positions
            tmpEvents = resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.pMat;
            for ithEvent=1:numel(tmpEvents)
                avgDecodePos = avgDecodePos + sum(tmpEvents{ithEvent}{1}.pMat, 2)';
                % figure; plot(sum(tmpEvents{ithEvent}{1}.pMat, 2))
            end
            % figure; plot(avgDecodePos)
        end
        
        
        if plotNetStruct
            figure; histogram(sum(network.cluster_mat, 1)); xlabel('n Clusters'); ylabel('Neurons (count)')
            histcounts(sum(network.cluster_mat, 1))
            % figure; histogram(sum(network.cluster_mat, 2)); xlabel('n Neurons'); ylabel('Clusters (count)')
        end
        
        % PF: plot and calculate score
        if calcScore
            % network = struct; linfields = {}; network.E_indices = E_indices; network.all_indices = 1:parameters.n;
            day = 1; epoch = 1; tetrode = 1; tr = 1;
            linfields{day}{epoch}{tetrode}{1}{tr}(:,1) = pfsim.gridxvals*100; % convert to cm
            for ithCell = network.E_indices
                linfields{day}{epoch}{tetrode}{ithCell}{tr}(:,5) = PFmat(ithCell,:);
            end
            PFscore = calculate_linfieldsScore(linfields, pfsim, pfsim, network, 'plotPFs', false);
            PFscores_all = [PFscores_all, PFscore];
        end
        
        
        if nEventsToPlot>0
            if isempty(resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.pvalue); continue; end
            
            % rvals_preplay = resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.rsquare(:,1); % pvalue
            % fitPvals_preplay = resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.FitpVal(:,1); % pvalue
            
            % PF Peak sequence
            PFmat_E = PFmat(E_indices,:);
            [peakRate,peakRateLocation] = max(PFmat_E, [], 2);

            if useMeanPFDensity % Sort by location of mean PF density
                PFexpectedLocation = sum( PFmat_E./sum(PFmat_E, 2) .*([1:size(PFmat_E, 2)]), 2); % units of space bins
                %PFexpectedLocation(peakRate<minPeakRate)=nan;
                [B,sortedCellIndsbyExpectedLocation] = sort(PFexpectedLocation, 'descend');
                PFpeaksSequence = sortedCellIndsbyExpectedLocation;
            else % Sort by location of PF peak
                %peakRateLocation(peakRate<minPeakRate)=nan;
                [B,sortedCellIndsbyPeakRateLocation] = sort(peakRateLocation, 'descend');
                PFpeaksSequence = sortedCellIndsbyPeakRateLocation;
            end

            eventLengths = resultsStruct(ithParam1, ithParam2, ithNet).results.eventLength;
            mat = resultsStruct(ithParam1, ithParam2, ithNet).results.ranksVec;
            x = mat./max(mat);
            x = x(:, eventLengths>=0.05); % temp line, since decoding has different min length
            x(peakRate<minPeakRate,:)=nan; % Exclude non high-rate cells

            [pvals_preplay_sorted, pvals_sorted ]= sort(pvals_preplay);

            % Cluster-wise coloring of raster
            % seqMat = x(PFpeaksSequence,pvals_sorted(1));
            clustMat = network.cluster_mat(:,E_indices)'; 
            clustColor = clustMat*[1:netParams.clusters]' ./ sum(clustMat, 2);

            % Accumulate best events
            for ithEvent = 1:size(pvals_preplay)
                if any(pvals_preplay(ithEvent) < bestEventpVals )
                    [~, minInd] = max(bestEventpVals);
                    bestEventpVals(minInd) = pvals_preplay(ithEvent);
                    bestEventjdVals(minInd) = maxJumps_preplay(ithEvent);
                    bestEventAbsrVals(minInd) =  sqrt(rsqrs_preplay(ithEvent));
                    bestEvents_decode{minInd} = resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.pMat{ithEvent}{1}.pMat;
                    bestEvents_relRank{minInd} = x(PFpeaksSequence,ithEvent);
                    bestEvents_relRankColor{minInd} = clustColor(PFpeaksSequence);
                end
            end
            nEvents = nEvents+numel(pvals_preplay);
            
            % Plot ithNet's best event
            if plotNetsBest
                %{
                %figure; scatter(1:parameters.n_E, x(PFpeaksSequence,pvals_sorted(1))); title(['Best event of network ', num2str(ithNet)])
                figure; imagesc(resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.pMat{pvals_sorted(1)}{1}.pMat )
                  xlabel('Time bin (10 ms)'); ylabel('Space bin (2 cm)'); title(['Best event of network ', num2str(ithNet)])
                figure; scatter(1:parameters.n_E, x(PFpeaksSequence,pvals_sorted(1)), [], clustColor(PFpeaksSequence), 'filled'); colorbar; 
                  xlabel('Neuron (sort by PF expected location)'); ylabel('Event relative rank'); title(['Best event of network ', num2str(ithNet)])
                %}
                  
                figure; scatter(x(PFpeaksSequence,pvals_sorted(1))', 1:parameters.n_E, 'k', '|', 'LineWidth', 5 )
                ylabel('Place cell (sorted)'); xlabel('Event relative rank'); title(['Best event of network ', num2str(ithNet)]); ylim([0, parameters.n_E])
                
                eventPmat = resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.pMat{pvals_sorted(1)}{1}.pMat;
                [yt, xt] = size(eventPmat);
                figure; imagesc([1:xt]*(tBinSz), [1:yt]*(pfsim.spatialBin*100), eventPmat)
                xlabel('Time (ms)'); ylabel('Position (cm)'); title(['Best event of network ', num2str(ithNet)])
                colormap(hot); title ''; colorbar off; caxis(([0, 0.25]))
                     
                  keyboard
            end
            
        end
        

    end

    
    %% Plot p-value matrix
    
    % myPlotSettings(3, 1.5) % for poster
    if plotPvalMat
        rvalThresh_vec = 0.0:0.1:0.9;
        jumpThres_vec =  0.1:0.1:1.0;

        op = zeros(numel(jumpThres_vec), numel(rvalThresh_vec));
        for ij = 1:numel(jumpThres_vec)
            allEventRs_pvalmat = sqrt(allEventRs);
            allShuffleRs_pvalmat = cellfun(@sqrt, allShuffleRs, 'UniformOutput', false);
            %allEventRs_pvalmat = allEventRs;
            %allShuffleRs_pvalmat = allShuffleRs;
            
            for ir = 1:numel(rvalThresh_vec)
                nActPass = mean( [allEventRs_pvalmat>rvalThresh_vec(ir)] & [allEventMaxJumps<jumpThres_vec(ij)]);
                nShuffPass = zeros(1, size(allShuffleMaxJumps{1}, 2));
                for ithShuf = 1:size(allShuffleMaxJumps{1}, 2)
                    nShuff_jump = cellfun(@(x) [x(ithShuf)<jumpThres_vec(ij)], allShuffleMaxJumps);
                    nShuff_rval = cellfun(@(x) [x(ithShuf)>rvalThresh_vec(ir)], allShuffleRs_pvalmat);
                    nShuffPass(ithShuf) = mean( nShuff_jump & nShuff_rval );
                end

                op(ij, ir) = 1 - mean(nActPass>nShuffPass);
                if [sum(nShuffPass)==0] & [sum(nActPass)==0]
                    op(ij, ir) = nan;
                end
            end
            ij
        end

        % myPlotSettings(4, 3, 2, 14, [], [], 2) % ppt format
        
        % Plot with linear colormap
        figure; imagesc(jumpThres_vec, rvalThresh_vec, op', 'AlphaData', ~isnan(op')); colorbar
        xlabel('<|Max Jump Distance|')
        ylabel('>|Correlation|')
        caxis([0, 0.1])
        
        
        % Plot with better colormap
        figure; 
        imagesc(jumpThres_vec, rvalThresh_vec, log10(op'), 'AlphaData', ~isnan(op'))
        % set(gca,'YDir','normal')
        cb = colorbar(); %cb.Label.String = cbLabel2;
        xlabel('<|Max Jump Distance|')
        ylabel('>|Correlation|')
        
        %title(analysisTitle)
        N = 256; n = N/2;
        cm = NaN(N,3);
        cm(:,1) = [ones(n,1);linspace(1,0,N-n)';];
        cm(:,2) = [linspace(0,1,n)';linspace(1,0,N-n)']; 
        cm(:,3) = [linspace(0,1,n)';ones(N-n,1)]; 
        set(gca,'clim',[log10(.05)*2 0])
        set(gcf,'colormap',cm)
        colorbar
        colorbar('Direction','reverse','Ticks',[log10(.005),log10(.05),log10(.5)],'TickLabels',[.005,.05,.5])
        hold on; 
        [xnan, ynan] = find(isnan(op));
        scatter(rvalThresh_vec(xnan), jumpThres_vec(ynan), 300, 'x')
    end

    %% Plot mean decode position
    figure; plot(avgDecodePos./sum(avgDecodePos))
    xlabel('Position (2 cm bins)'); ylabel('Probability Density')
    
    
    %% Plot ECDF of r values
    
    % myPlotSettings(3, 1.5) % for poster
    
    myPlotSettings(4, 2.75, 1.5, 14, [], [], 1.25) % ppt format
    
    allEventRs_ecdf = (allEventRs);
    rvals_shuff_ecdf= (vertcat(allShuffleRs{:}));
                
    figure; hold on; ecdf(sqrt(allEventRs_ecdf)); ecdf(sqrt(rvals_shuff_ecdf)); %ecdf(sqrt(rvals_shuff_ecdf)); 
    % figure; hold on; ecdf(allEventRs_ecdf); ecdf(rvals_ecdf)
    legend({'Preplays', 'Shuffles'}, 'Location', 'Best')
    xlabel('|Correlation|'); ylabel('Cumulative proportion');
    [H,P,KSSTAT] = kstest2(sqrt(allEventRs_ecdf), sqrt(rvals_shuff_ecdf) )
    title([variedParam(1).name, '=', num2str(variedParam(1).range(ithParam1)), ' ', ...
        variedParam(2).name, '=', num2str(variedParam(2).range(ithParam2)), ...
        ' pval=', num2str(P), ...
        ' nEvents=', num2str(numel(allEventRs_ecdf))])

    title ''
    legend({'Preplay', 'Shuffle'}, 'Location', 'Best'); legend boxoff
    
    
    h = get(gca,'children');
%    set(h(1), 'LineWidth', 1, 'color', [0.7, 0, 0], 'LineStyle', '-.')
    set(h(1), 'LineWidth', 1, 'color', [1, 0.3, 0.3], 'LineStyle', '-.')

    set(h(2), 'LineWidth', 1, 'color', [0, 0, 0.6])
    %set(h(2), 'LineWidth', 1, 'color', 'k', 'LineStyle', '--'); set(h(3), 'LineWidth', 1, 'color', [0, 0, 0.6])


    %% Plot best sequences
    
    % Decodes
    myPlotSettings(7, 6)
    % myPlotSettings(6, 4) % for poster, main result
    % myPlotSettings(3.5, 5.5) % for poster, secondary results
     myPlotSettings(4.5, 3, 2, 12, [], [], 2) % ppt format
    
    if nEventsToPlot>0
        [pvals_sorted, eventSortInd] = sort(bestEventpVals);
        figure; tfig = tiledlayout('flow', 'Padding', 'none', 'TileSpacing', 'compact');
        title(tfig, ['Best ', num2str(numel(bestEvents_decode)), ' of ', num2str(nEvents), ' events'])
        for ithEventToPlot = eventSortInd'
            nexttile
            
            [yt, xt] = size(bestEvents_decode{ithEventToPlot});
            imagesc([1:xt]*(tBinSz), [1:yt]*(pfsim.spatialBin*100), bestEvents_decode{ithEventToPlot})
            colormap hot
            title(['pval=', num2str(bestEventpVals(ithEventToPlot))], 'FontWeight', 'normal')
            title(['r=', num2str(bestEventAbsrVals(ithEventToPlot), 2), '; jd=', num2str(bestEventjdVals(ithEventToPlot), 2)], 'FontWeight', 'normal')            
            % caxis(([0, 0.75*max(bestEvents_decode{ithEventToPlot}, [], 'all')]))
            caxis(([0, 0.25]))
             
            set(gca,'ytick',[])
            %yt = yticks; yticklabels(yt*(pfsim.spatialBin*100)); ylabel('Position (cm)')
            %xt = xticks; xticklabels(xt*(tBinSz)); xlabel('Time (ms)')
            %xticks('auto')
        end
    end
    title(tfig, '')
    
    
    myPlotSettings
    
    % Relative ranks (corresponding to decodes)
    if nEventsToPlot>0
        [pvals_sorted, eventSortInd] = sort(bestEventpVals);
        figure; tfig = tiledlayout('flow');
        title(tfig, ['Best ', num2str(numel(bestEvents_relRank)), ' of ', num2str(nEvents), ' events'])
        for ithEventToPlot = eventSortInd'
            nexttile
            eventSeq = bestEvents_relRank{ithEventToPlot};
            scatter(1:sum(~isnan(eventSeq)), eventSeq(~isnan(eventSeq)), [], bestEvents_relRankColor{ithEventToPlot}(~isnan(eventSeq)), 'filled')
            %scatter(1:parameters.n_E, bestEvents_relRank{ithEventToPlot})
            title(['pval=', num2str(bestEventpVals(ithEventToPlot))])
            %caxis(([0, 0.5*max(bestEvents{ithEventToPlot}, [], 'all')]))
        end
    end
    

    %% Plot combined place fields
    
    % myPlotSettings(4, 3, 2, 14, [], [], 2) % ppt format
    
    row_all_zeros1 = find(all( PFmatE_all==0, 2)) ;
    row_n_all_zeros1 = find(~all( PFmatE_all==0, 2)) ;
    [peakRate,peakRateLocation] = max(squeeze(PFmatE_all(row_n_all_zeros1,:)), [], 2);
    [B,sortedCellIndsbyPeakRateLocation] = sort(peakRateLocation, 'descend');
    PFpeaksSequence = [row_n_all_zeros1(sortedCellIndsbyPeakRateLocation); row_all_zeros1];
    [peakRate, peakRateLocation_all] = max(PFmatE_all, [], 2);

    normRates = 1;
    if normRates
        rateDenomAll = max(PFmatE_all(PFpeaksSequence,:), [], 2);
        caxmax = 1;
    else
        rateDenomAll = 1;
        caxmax = max(PFmatE_all, [], 'all');
    end
    
    PFdatatoplot = PFmatE_all(PFpeaksSequence(rateDenomAll>minPeakRate),:)./rateDenomAll(rateDenomAll>minPeakRate);
    figure; imagesc( fliplr(PFdatatoplot(:,5:end)) ); title('All nets'); colorbar; caxis([0, caxmax])
    xlabel('Position (2 cm bin)'); ylabel('Cell (sorted)');
    
    colormap(jet); title ''; colorbar off
    xt = xticks; xticklabels(xt*(pfsim.spatialBin*100)); xlabel('Position (cm)')
    % set(gca,'xtick',[]); xlabel('')
    
    %% Plot cluster-wise place fields
    if plotClusterPFs
        figure; tiledlayout(parameters.clusters,1);
        singularMembership = 1;
        for ithCluster = parameters.clusters:-1:1
            if singularMembership
                isClusterMember = [sum(cluster_matAll(PFpeaksSequence(rateDenomAll>minPeakRate), :), 2)==1] .* [cluster_matAll(PFpeaksSequence(rateDenomAll>minPeakRate),ithCluster)==1];
            else
                isClusterMember = [cluster_matAll(PFpeaksSequence(rateDenomAll>minPeakRate),ithCluster)==1];
            end
            clusterPF = isClusterMember .* PFmatE_all(PFpeaksSequence(rateDenomAll>minPeakRate),:)./rateDenomAll(rateDenomAll>minPeakRate);
            zeroInds = all(clusterPF==0, 2);
            nexttile; imagesc( clusterPF(~zeroInds,:) ); 
            %xaxis('off')
            set(gca,'xtick',[])

            % title('All nets'); colorbar; caxis([0, caxmax])
            % xlabel('Position (2 cm bin)'); ylabel('Cell (sorted)');
        end
    end
    
    
    %% Extra plots      
    if plotExtraPlots
        % Example network PFs    
        for ithNet = 1%:10
            cellsToPlot = [rateDenomAll>minPeakRate] & [netInd_all(PFpeaksSequence)==ithNet];
            figure; imagesc( PFmatE_all(PFpeaksSequence(cellsToPlot),:)./rateDenomAll(cellsToPlot)); title('Example net'); colorbar; caxis([0, caxmax])
            xlabel('Position (2 cm bin)'); ylabel('Cell (sorted)');
        end

        % 'Peak Rate'
        % myPlotSettings(3, 1.5) % for poster
        peakRate= max(PFmatE_all, [], 2);
        figure; histogram(peakRate(peakRate>minPeakRate)); xlabel('Place field peak (Hz)'); ylabel('Place cells (count)');

        % 'kstest'
        ksstat_vec = [];
        for i = 1:size(PFmatE_all, 1)
            [~,p,ksstat,~] = kstest( ( PFmatE_all(i,:)-mean(PFmatE_all(i,:), 2) )./(std(PFmatE_all(i,:), [], 2)+eps  ) );
            ksstat_vec(i) = ksstat;
        end   
        figure; histogram(ksstat_vec(peakRate>minPeakRate)); xlabel('kstest stat'); ylabel('Place cells (count)');

        % 'sparsity'                
        cellSparsity =  mean( PFmatE_all>[0.25*max(PFmatE_all, [], 2)], 2 );
        figure; histogram(cellSparsity(peakRate>minPeakRate)); xlabel('PF sparsity'); ylabel('Place cells (count)');

        % 'information'
        spatialInfo = nanmean( [PFmatE_all./mean(PFmatE_all, 2)] .* log(( PFmatE_all+eps )./mean(PFmatE_all, 2) ), 2 );
        figure; histogram(spatialInfo(peakRate>minPeakRate)); xlabel('Spatial information (bits/s)'); ylabel('Place cells (count)');

        
        figure; scatter(cellSparsity, spatialInfo); xlabel('PF sparsity'); ylabel('Spatial info.')
        figure; scatter(peakRate, spatialInfo); xlabel('Peak rate (Hz)'); ylabel('Spatial info.')
        
        mdl1 = fitlm(nClustMemb_all, spatialInfo); figure; plot(mdl1); xlabel('n cluster membership'); ylabel('Spatial information'); 
        title(['pval=', num2str(mdl1.Coefficients.pValue(2), 2), ', rsqr=', num2str(mdl1.Rsquared.Ordinary, 2)]); legend off; mdl1
        mdl2 = fitlm(nClustMemb_all, cellSparsity); figure; plot(mdl2); xlabel('n cluster membership'); ylabel('Spatial sparsity'); 
        title(['pval=', num2str(mdl2.Coefficients.pValue(2), 2), ', rsqr=', num2str(mdl2.Rsquared.Ordinary, 2)]); legend off; mdl2
        mdl3 = fitlm(nClustMemb_all, peakRate); figure; plot(mdl3); xlabel('n cluster membership'); ylabel('Peak rate (Hz)'); 
        title(['pval=', num2str(mdl3.Coefficients.pValue(2), 2), ', rsqr=', num2str(mdl3.Rsquared.Ordinary, 2)]); legend off; mdl3

        cellsToPlot = [spatialInfo>1.87]; % [rateDenomAll>2];
        figure; imagesc( PFmatE_all(PFpeaksSequence(cellsToPlot),:)./rateDenomAll(cellsToPlot)); title('All nets'); colorbar; caxis([0, caxmax])
        xlabel('Position (2 cm bin)'); ylabel('Cell (sorted)');
        
        
        % 'Example best place fields'
        [~, SIind] = sort( ksstat_vec' .* [peakRateLocation_all>10&peakRateLocation_all<40] .* [peakRate>4] .* [cellSparsity<0.5], 'descend' );
        % [~, SIind] = sort( spatialInfo' .* [peakRateLocation_all>10&peakRateLocation_all<40] .* [meanPeakRate>4] .* [cellSparsity>0.5], 'descend' );
        figure; plot(1:size(PFmatE_all, 2), PFmatE_all(SIind(1:4),:))
        xlabel('Location (2 cm bins)'); ylabel('Firing rate (Hz)'); title('Good example place fields')
        figure; plot(1:size(PFmatE_all, 2), PFmatE_all(SIind(5:8),:))
        xlabel('Location (2 cm bins)'); ylabel('Firing rate (Hz)'); title('Good example place fields')
        figure; plot(1:size(PFmatE_all, 2), PFmatE_all(SIind(9:12),:))
        xlabel('Location (2 cm bins)'); ylabel('Firing rate (Hz)'); title('Good example place fields')
    end
    
    
    %% Plot PF score, if it was calculated
    if calcScore
        figure; histogram(PFscores_all); xlabel('Place field score (lower is better)'); ylabel('Network (count)');
    end
    
    % Plots for comparison to Shin et al., 2019
    if plotPFStats
        % 'Peak Rate'
        peakRate= max(PFmatE_all, [], 2);
        figure; histogram(peakRate(peakRate>minPeakRate), 10, 'Normalization', 'probability'); xlabel('Place field peak (Hz)'); ylabel('Place cells (Prob.)');
        xlim([-2, 45]); ylim([0, 0.4]); box off
        % 'sparsity'                
        cellSparsity =  mean( PFmatE_all>[0.25*max(PFmatE_all, [], 2)], 2 );
        figure; histogram(cellSparsity(peakRate>minPeakRate), 10, 'Normalization', 'probability'); xlabel('Spatial sparsity'); ylabel('E cells (count)'); ylabel('Place cells (Prob.)');
        xlim([-0.05, 1.05]); ylim([0, 0.3]); box off
        % 'information'
        spatialInfo = nanmean( [PFmatE_all./mean(PFmatE_all, 2)] .* log(( PFmatE_all+eps )./mean(PFmatE_all, 2) ), 2 );
        figure; histogram(spatialInfo(peakRate>minPeakRate), 10, 'Normalization', 'probability'); xlabel('Spatial Information (bits)'); ylabel('E cells (count)'); ylabel('Place cells (Prob.)');
        xlim([-0.05, 3.45]); ylim([0, 0.3]); box off
        
        % myPlotSettings(4, 2, 2, 14, [], [], 2) % ppt format
        xx = 0:0.1:10;
        yy = [gaussmf(xx,[1 2]); gaussmf(xx,[1 3]); gaussmf(xx,[1 4]); gaussmf(xx,[1 5]); gaussmf(xx,[1 6]); gaussmf(xx,[1 7]); gaussmf(xx,[1 8]);]
        yy = [gaussmf(xx,[1 3]); gaussmf(xx,[1 4]); gaussmf(xx,[1 5]);]
        figure; plot(xx,yy)
        xlabel('Position'); ylabel('Input rate (Hz)')
        box off; set(gca,'xtick',[]); set(gca,'ytick',[])
        
        xx = 0:0.5:1;
        yy = [ ([6, 6, 6].*xx(end)); (6.*xx); (6.*fliplr(xx))];
        figure; plot(xx,yy')
        xlabel('Position'); ylabel('Input rate (kHz)')
        box off;     
        % set(gca,'xtick',[]); set(gca,'ytick',[]); ylabel('Input rate (Hz)')
    end
    
    %% Plot new analyses for Pizza talk questions/supplement
    
    % Correlation of event length with decode correlation
    mdl_preplay = fitlm(allEventLengths, allEventRs); figure; plot(mdl_preplay); xlabel('Event time bins (count)'); ylabel('Decode correlation (|r|)'); 
    title(['pval=', num2str(mdl_preplay.Coefficients.pValue(2), 2), ', rsqr=-', num2str(sqrt(mdl_preplay.Rsquared.Ordinary), 2)]); legend off; mdl_preplay
    figure; scatterhist(allEventLengths, allEventRs); xlabel('Event time bins (count)'); ylabel('Decode correlation (|r|)'); 
    
    mdl_shuff = fitlm(allShuffleLengths, vertcat(allShuffleRs{:}) ); figure; plot(mdl_shuff); xlabel('Event time bins (count)'); ylabel('Decode correlation (|r|)'); 
    title(['pval=', num2str(mdl_shuff.Coefficients.pValue(2), 2), ', rsqr=-', num2str(sqrt(mdl_shuff.Rsquared.Ordinary), 2)]); legend off; mdl_shuff
    figure; scatterhist(allShuffleLengths, vertcat(allShuffleRs{:})); xlabel('Event time bins (count)'); ylabel('Decode correlation (|r|)'); 
    
    % 
    
    
end

runTime = toc;
disp([ 'Runtime: ', datestr(datenum(0,0,0,0,0,runTime),'HH:MM:SS') ])

