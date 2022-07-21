
%% Plot combined place fields at a particular parameter point
% load('C:\Users\Jordan\Box\Data\Replay project\RandNet\temp, new PF sim code grids\results_2022-07-12T06-58.mat')
% load('C:\Users\Jordan\Box\Data\Replay project\RandNet\temp, new PF sim code grids\results_2022-07-21T15-05.mat')

if isfolder('functions')
    addpath('functions')
else
    addpath( fullfile( '..', 'functions' ) )
end

% paramSetInds = [3, 3; 3 4; 4,2]
% paramSetInds = combvec([1:size(resultsStruct, 1)], [1:size(resultsStruct, 2)])'
paramSetInds = combvec([2], [2])'

parameters.n
parameters.del_G_sra
variedParam(1).name
variedParam(2).name


minPeakRate = 2; % minimum peak PF rate to be considered a place cell
calcScore = 0
plotExtraPlots = 0 % if 1, plot place fields of every network
nEventsToPlot = 25
useMeanPFDensity = true
plotNetsBest = false

%% Loop over parameter sets
rng(1); tic
for ithParamSet = 1:size(paramSetInds, 1)
    
    ithParamSet
	ithParam1 = paramSetInds(ithParamSet,1);
	ithParam2 = paramSetInds(ithParamSet,2);
    
        
    temp = nan(1, num_nets);
    
    PFmatE_all = [];
    netInd_all = []; % Index of which net each cell came from
    nClustMemb_all = []; % number of clusters each cell is a member of
    PFscores_all = [];
    bestEventpVals = inf(nEventsToPlot, 1); % nEventsToPlot best event pvals
    bestEvents = cell(nEventsToPlot, 1);
    nEvents = 0;
    cluster_matAll = [];
    for ithNet = 1:size(resultsStruct, 3)

        % Get matrix of PFs
        PFmat = PFresultsStruct(ithParam1, ithParam2, ithNet).results{1}.linfields;
        E_indices = PFresultsStruct(ithParam1, ithParam2, ithNet).results{1}.E_indices;
        netSeed = PFresultsStruct(ithParam1, ithParam2, ithNet).results{1}.netSeed

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
            
            % Accumulate best events
            pvals_preplay = resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.pvalue(:,1); % pvalue
            % rvals_preplay = resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.rsquare(:,1); % pvalue
            % fitPvals_preplay = resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.FitpVal(:,1); % pvalue
            for ithEvent = 1:size(pvals_preplay)
                if any(pvals_preplay(ithEvent) < bestEventpVals )
                   [~, minInd] = max(bestEventpVals);
                   bestEventpVals(minInd) = pvals_preplay(ithEvent);
                   bestEvents{minInd} = resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.pMat{ithEvent}{1}.pMat;
                end
            end
            nEvents = nEvents+numel(pvals_preplay);
            
            % Plot ithNet's best event
            if plotNetsBest
                % PF Peak sequence
                PFmat_E = PFmat(E_indices,:);
                [peakRate,peakRateLocation] = max(PFmat_E, [], 2);

                if 1 useMeanPFDensity % Sort by location of mean PF density
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
                x = x(:, eventLengths>0.05); % temp line, since decoding has different min length
                x(peakRate<minPeakRate,:)=nan; % Exclude non high-rate cells

                [pvals_preplay_sorted, pvals_sorted ]= sort(pvals_preplay);
                figure; scatter(1:parameters.n_E, x(PFpeaksSequence,pvals_sorted(1))); title(['Best event of network ', num2str(ithNet)])
                figure; imagesc(resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.pMat{pvals_sorted(1)}{1}.pMat )
            end
            
        end
        

    end

    
    %% Plot best sequences

    if nEventsToPlot>0
        [pvals_sorted, eventSortInd] = sort(bestEventpVals);
        figure; tfig = tiledlayout('flow');
        title(tfig, ['Best ', num2str(numel(bestEvents)), ' of ', num2str(nEvents), ' events'])
        for ithEventToPlot = eventSortInd'
            nexttile
            imagesc(bestEvents{ithEventToPlot})
            title(['pval=', num2str(bestEventpVals(ithEventToPlot))])
            %caxis(([0, 0.5*max(bestEvents{ithEventToPlot}, [], 'all')]))
        end
    end

    keyboard
    

    %% Plot combined place fields
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
    
    figure; imagesc( PFmatE_all(PFpeaksSequence(rateDenomAll>minPeakRate),:)./rateDenomAll(rateDenomAll>minPeakRate)); title('All nets'); colorbar; caxis([0, caxmax])
    xlabel('Position (2 cm bin)'); ylabel('Cell (sorted)');
    % colormap(flipud(bone))
    
    %% Plot cluster-wise place fields
    figure; tiledlayout(parameters.clusters,1);
    singularMembership = 0;
    for ithCluster = parameters.clusters:-1:1
        if singularMembership
            isClusterMember = [sum(cluster_matAll(PFpeaksSequence(rateDenomAll>minPeakRate), :), 2)==1] .* [cluster_matAll(PFpeaksSequence(rateDenomAll>minPeakRate),ithCluster)==1];
        else
            isClusterMember = [cluster_matAll(PFpeaksSequence(rateDenomAll>minPeakRate),ithCluster)==1];
        end
        clusterPF = isClusterMember .* PFmatE_all(PFpeaksSequence(rateDenomAll>minPeakRate),:)./rateDenomAll(rateDenomAll>minPeakRate);
        zeroInds = all(clusterPF==0, 2);
        nexttile; imagesc( clusterPF(~zeroInds,:) ); %title('All nets'); colorbar; caxis([0, caxmax])
        %xlabel('Position (2 cm bin)'); ylabel('Cell (sorted)');
    end

    %% Extra plots      
    if plotExtraPlots
        % Example network PFs    
        cellsToPlot = [rateDenomAll>minPeakRate] & [netInd_all(PFpeaksSequence)==10];
        figure; imagesc( PFmatE_all(PFpeaksSequence(cellsToPlot),:)./rateDenomAll(cellsToPlot)); title('Example net'); colorbar; caxis([0, caxmax])
        xlabel('Position (2 cm bin)'); ylabel('Cell (sorted)');

        % 'Peak Rate'
        meanPeakRate= max(PFmatE_all, [], 2);
        figure; histogram(meanPeakRate(rateDenomAll>minPeakRate)); xlabel('Peak PF rate (Hz)'); ylabel('E cells (count)');

        % 'kstest'
        ksstat_vec = [];
        for i = 1:size(PFmatE_all, 1)
            [~,p,ksstat,~] = kstest( ( PFmatE_all(i,:)-mean(PFmatE_all(i,:), 2) )./(std(PFmatE_all(i,:), [], 2)+eps  ) );
            ksstat_vec(i) = ksstat;
        end   
        figure; histogram(ksstat_vec(rateDenomAll>minPeakRate)); xlabel('kstest stat'); ylabel('E cells (count)');

        % 'sparsity'                
        cellSparsity = mean( PFmatE_all<=[0.25*max(PFmatE_all, [], 2)], 2 );
        figure; histogram(cellSparsity(rateDenomAll>minPeakRate)); xlabel('PF sparsity'); ylabel('E cells (count)');

        % 'information'
        spatialInfo = nanmean( [PFmatE_all./mean(PFmatE_all, 2)] .* log(( PFmatE_all+eps )./mean(PFmatE_all, 2) ), 2 );
        figure; histogram(spatialInfo(rateDenomAll>minPeakRate)); xlabel('PF information'); ylabel('E cells (count)');

        %{
        figure; scatter(cellSparsity, spatialInfo)
        figure; scatter(meanPeakRate, spatialInfo)
        
        mdl1 = fitlm(nClustMemb_all, spatialInfo); figure; plotAdded(mdl1); xlabel('n cluster membership'); ylabel('Spatial info.'); mdl1
        mdl2 = fitlm(nClustMemb_all, cellSparsity); figure; plot(mdl2); xlabel('n cluster membership'); ylabel('PF sparsity'); mdl2
        mdl3 = fitlm(nClustMemb_all, meanPeakRate); figure; plot(mdl3); xlabel('n cluster membership'); ylabel('Peak rate (Hz)'); mdl3
        

        cellsToPlot = [spatialInfo>1.87]; % [rateDenomAll>2];
        figure; imagesc( PFmatE_all(PFpeaksSequence(cellsToPlot),:)./rateDenomAll(cellsToPlot)); title('All nets'); colorbar; caxis([0, caxmax])
        xlabel('Position (2 cm bin)'); ylabel('Cell (sorted)');
        %}
        
        % 'Example best place fields'
        [~, SIind] = sort( ksstat_vec' .* [peakRateLocation_all>10&peakRateLocation_all<40] .* [meanPeakRate>4] .* [cellSparsity>0.5], 'descend' );
        % [~, SIind] = sort( spatialInfo' .* [peakRateLocation_all>10&peakRateLocation_all<40] .* [meanPeakRate>4] .* [cellSparsity>0.5], 'descend' );
        figure; plot(1:size(PFmatE_all, 2), PFmatE_all(SIind(1:4),:))
        xlabel('Location (2 cm bins)'); ylabel('Firing rate (Hz)'); title('Example place fields')
        
    end
    
    %% Plot PF score, if it was calculated
    if calcScore
        figure; histogram(PFscores_all); xlabel('Place field score (lower is better)'); ylabel('Network (count)');
    end
    
end

runTime = toc;
disp([ 'Runtime: ', datestr(datenum(0,0,0,0,0,runTime),'HH:MM:SS') ])

