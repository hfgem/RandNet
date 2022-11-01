
if isfolder('functions')
    addpath('functions')
else
    addpath( fullfile( '..', 'functions' ) )
end

if ispc
    addpath('C:\Users\Jordan\Box\Data\RandNet-Data\temp, new PF sim code grids')
elseif ismac
    addpath(['/Users/jordan/Library/CloudStorage/Box-Box/Data/RandNet-Data/temp, new PF sim code grids'])
else
    disp('error')
end


 load('results_2022-07-22T18-29.mat') % Primary clustersXmnc grid (smaller mnc values)
% load('results_2022-07-13T11-44.mat') % clustersXmnc grid up to (5,5)


%% Plot preplay decoding results across grid search
%
% resultsStruct(param1Ind, param2Ind, ithNet)
% PFresultsStruct(param1Ind, param2Ind, ithNet)
%
% num_nets 
% variedParam
%

%{
load('results\randnet_PF_param_tests\results_2022-06-29T14-17.mat')
%}

parameters.n
parameters.del_G_sra
variedParam(1).name
variedParam(2).name

combineNetData = 1 
plotAllCDFs = 0
plotCombinedCDFs = 0

maxNEvents = inf % downsample number of replay events for kstest to maxNEvents, use inf for no downsampling
% maxNEvents = 500 

useWeightedDecode = 0; % slope and R2 for correlations by either peak prob or weighted prob

%downSample = 0
%downSampleFracEvents = 0.5;

xParamvec = variedParam(1).range;
xName = variedParam(1).name;
yParamvec = variedParam(2).range;
yName = variedParam(2).name;

if useWeightedDecode
    analysisTitle = 'KS-test, weighted decode corrs';
else
    analysisTitle = 'KS-test, peak pos. decode corrs';
end
if combineNetData
    cbLabel1 = 'Combined p-val';
    cbLabel2 = 'Combined KS-stat';
else
    cbLabel1 = 'Median p-val';
    cbLabel2 = 'Median KS-stat';
end

op = nan(3, numel(xParamvec), numel(yParamvec));

allNetPvals = zeros(numel(xParamvec), numel(yParamvec), num_nets);
allNetKSstat = zeros(numel(xParamvec), numel(yParamvec), num_nets);
allNetMedianDiff = zeros(numel(xParamvec), numel(yParamvec), num_nets);
allNetGroupID = nan(numel(xParamvec), numel(yParamvec), num_nets);

allNetDecodeSlopes = zeros(numel(xParamvec), numel(yParamvec), num_nets); % Mean decoded slope of each network for each parameter point
allNetDecodeSlopes_zscored = zeros(numel(xParamvec), numel(yParamvec), num_nets); % Mean decoded slope of each network for each parameter point

allNetDecodeEntropy = zeros(numel(xParamvec), numel(yParamvec), num_nets); % Mean decoded slope of each network for each parameter point

allNetDecodeVar = zeros(numel(xParamvec), numel(yParamvec), num_nets); % Mean decoded slope of each network for each parameter point

figure(515); hold on
imagesc(xParamvec, yParamvec, squeeze(op(1,:,:))', 'AlphaData', ~isnan(squeeze(op(1,:,:))')); drawnow
set(gca,'YDir','normal')
cb = colorbar();
xlabel(xName,'Interpreter','none')
ylabel(yName,'Interpreter','none')

rng(2)
% gcp
tic
for ithParam1 = 1:size(resultsStruct, 1)
    for ithParam2 = 1:size(resultsStruct, 2)
        
        temp = nan(1, num_nets);
        allRvecs = [];
        allRvecs_shuff = [];
        nEvents_accum = 0;
        nSigEvents_accum = 0;
        for ithNet = 1:size(resultsStruct, 3)
            
            % if ~isempty(resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.pMat)
            if [~isempty(resultsStruct(ithParam1, ithParam2, ithNet).results) && ...
                isfield(resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory, 'pMat') && ...
                ~isempty(resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.pMat)]

                % keyboard
                %pvals = resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.pvalue(:,1);
                %figure; histogram(pvals, 10)

                if useWeightedDecode
                    allshuff_rvals = vertcat(resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.shuffle_weightedR2{:});
                    allshuff_rvals = allshuff_rvals(:,1); % take just forward traj
                    rvals_shuffle = vertcat([allshuff_rvals{:,1}]');
                    
                    rvals_preplay = resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.weightedR2(:,1);
                
                    pvals_preplay = nan(size(rvals_preplay));
                    for ithEvent=1:numel(rvals_preplay)
                        pvals_preplay(ithEvent) = mean(rvals_preplay(ithEvent)<allshuff_rvals{ithEvent});
                    end
                    
                    slopes_preplay = resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.weightedSlope(:,1);
                    slopes_shuffle = vertcat(resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.shuffle_weightedSlope{:});
                    slopes_shuffle = slopes_shuffle(:,1); % take just forward traj
                    slopes_shuffle = vertcat([slopes_shuffle{:,1}]');
                else
                    allshuff_rvals = vertcat(resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.shuffle_rsquare{:});
                    allshuff_rvals = allshuff_rvals(:,1); % take just forward traj
                    rvals_shuffle = vertcat(allshuff_rvals{:,1});
                    
                    rvals_preplay = resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.rsquare(:,1);
                
                    % Equal to frac of shuffle events with greater R2 than actual event's R2
                    pvals_preplay = resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.pvalue(:,1);
                    
                    slopes_preplay = resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.slopes(:,1);
                    try % Struct field added to later simulation
                    slopes_shuffle = vertcat(resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.shuffle_slopes{:});
                    slopes_shuffle = slopes_shuffle(:,1); % take just forward traj
                    slopes_shuffle = cellfun(@transpose,slopes_shuffle,'UniformOutput',false); 
                    slopes_shuffle = vertcat([slopes_shuffle{:,1}]');
                    catch
                        slopes_shuffle = nan;
                    end
                end
                
                if plotAllCDFs
                    figure; hold on; ecdf(rvals_preplay); ecdf(rvals_shuffle); legend({'Preplays', 'Shuffles'}, 'Location', 'Best')
                end
                
                try
                    [H,P,KSSTAT] = kstest2(rvals_preplay, rvals_shuffle);
                    temp(1, ithNet) = P;
                    temp(2, ithNet) = KSSTAT;
                catch 
                    disp('Error, no events')
                    temp(1, ithNet) = nan;
                    temp(2, ithNet) = nan;
                end
                
                allNetPvals(ithParam1, ithParam2, ithNet) = P;
                allNetKSstat(ithParam1, ithParam2, ithNet) = KSSTAT;

                
                %[F_cdfAactual,X_cdfAactual] = ecdf(rvals_preplay); [F_cdfShuff,XcdfShuff] = ecdf(rvals_shuffle);
                %F_cdfAactualInterp = interp1(X_cdfAactual(2:end), F_cdfAactual(2:end), XcdfShuff(2:end));
                %allNetAUCdiff(ithParam1, ithParam2, ithNet) = trapz(F_cdfAactualInterp,XcdfShuff(2:end)) - trapz(F_cdfShuff,XcdfShuff);
                
                allNetMedianDiff(ithParam1, ithParam2, ithNet) = median(rvals_preplay) - median(rvals_shuffle);
                allNetGroupID(ithParam1, ithParam2, ithNet) = sub2ind( size(resultsStruct, [1 2]), ithParam1, ithParam2) ;
                
                % if allNetMedianDiff(ithParam1, ithParam2, ithNet)>0.04; keyboard; end
                
                allRvecs = [allRvecs, rvals_preplay'];
                allRvecs_shuff = [allRvecs_shuff, rvals_shuffle'];
                
                nEvents_accum = nEvents_accum + numel(pvals_preplay);
                nSigEvents_accum = nSigEvents_accum + sum(pvals_preplay<0.05);
                
                %%
                
                % slopes_shuffle = vertcat(resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.shuffle_slopes{:});
                % slopes_shuffle = slopes_shuffle(:,1); % take just forward traj
                allNetDecodeSlopes(ithParam1, ithParam2, ithNet) = mean(abs(slopes_preplay));
                allNetDecodeSlopes_zscored(ithParam1, ithParam2, ithNet) = [mean(abs(slopes_preplay)) - mean(abs(slopes_shuffle))] / [std(abs(slopes_shuffle))];

                entropy_preplay = resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.Entropy(:,1);
                allNetDecodeEntropy(ithParam1, ithParam2, ithNet) = mean(entropy_preplay);
                
                temp = 0;
                for ithEvent = 1:numel(resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.pMat)
                    temp = temp + mean(var(resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.pMat{ithEvent}{1}.pMat, [], 1));
                end
                allNetDecodeVar(ithParam1, ithParam2, ithNet) = temp./ numel(resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.pMat);
                                

                temp(3, ithNet) = numel(pvals_preplay) / sum(pvals_preplay<0.05) ;

            else
                temp(1, ithNet) = nan;
                temp(2, ithNet) = nan;
                temp(3, ithNet) = nan;
            end
                
        end
        
        if combineNetData
            if ~isempty(allRvecs)
                allRvecs = allRvecs(randperm(numel(allRvecs))); % permute, to mix networks' events
                allRvecs = allRvecs(1:min(numel(allRvecs), maxNEvents));
                numel(allRvecs)
                
                [~,p_kstest,KSSTAT] = kstest2(allRvecs, allRvecs_shuff);
                op(1, ithParam1, ithParam2) = p_kstest;
                op(2, ithParam1, ithParam2) = KSSTAT;
                op(3, ithParam1, ithParam2) = nSigEvents_accum/nEvents_accum;
            else
                p_kstest = nan;
                op(1, ithParam1, ithParam2) = nan;
                op(2, ithParam1, ithParam2) = nan;
                op(3, ithParam1, ithParam2) = nan;
            end
            
            if plotCombinedCDFs
                try
                figure; hold on; ecdf(allRvecs); ecdf(allRvecs_shuff); legend({'Preplays', 'Shuffles'}, 'Location', 'Best')
                title(['ithParam1=', num2str(ithParam1) ' ithParam2=', num2str(ithParam2), ' pval=', num2str(p_kstest) ])
                % [Fcdf,Xcdf] = ecdf(allRvecs); [Fcdfshuff,Xcdfshuff] = ecdf(allRvecs_shuff);
                N = ecdfhist(Fcdf, Xcdf);
                Nshuff = ecdfhist(Fcdfshuff, Xcdf);

                end
            end
            
        else
            op(1, ithParam1, ithParam2) = median(temp(1,:));
            %op(1, ithParam1, ithParam2) = min(temp(1,:));
            op(2, ithParam1, ithParam2) = nanmean(temp(2,:));
            op(3, ithParam1, ithParam2) = nanmean(temp(3,:));
        end
        

        %op(1, ithParam1, ithParam2) =  min(temp(1,:)); % nanmean(temp(1,:));
        %op(2, ithParam1, ithParam2) = nanmean(temp(2,:));
    end
    
    figure(515)
    imagesc(xParamvec, yParamvec, squeeze(op(1,:,:))', 'AlphaData', ~isnan(squeeze(op(1,:,:))')); drawnow
    
end
runTime = toc;
disp([ 'Runtime: ', datestr(datenum(0,0,0,0,0,runTime),'HH:MM:SS') ])
% disp( duration(0, 0, runTime) )


% myPlotSettings(3, 1.5) % for poster
% myPlotSettings(4, 3, 2, 14, [], [], 2) % ppt format


figure; 
imagesc(xParamvec, yParamvec, squeeze(op(1,:,:))', 'AlphaData', ~isnan(squeeze(op(1,:,:))'))
set(gca,'YDir','normal')
cb = colorbar(); cb.Label.String = cbLabel1;
xlabel(xName,'Interpreter','none')
ylabel(yName,'Interpreter','none')
title(analysisTitle)

figure; 
imagesc(xParamvec, yParamvec, squeeze(op(2,:,:))', 'AlphaData', ~isnan(squeeze(op(2,:,:))'))
set(gca,'YDir','normal')
cb = colorbar(); cb.Label.String = cbLabel2;
xlabel(xName,'Interpreter','none')
ylabel(yName,'Interpreter','none')
title(analysisTitle)


figure; 
imagesc(xParamvec, yParamvec, squeeze(op(3,:,:))', 'AlphaData', ~isnan(squeeze(op(3,:,:))'))
set(gca,'YDir','normal')
cb = colorbar(); cb.Label.String = 'Frac Sig';
xlabel(xName,'Interpreter','none')
ylabel(yName,'Interpreter','none')
title(analysisTitle)


% caxis([prctile(op, 2.5, 'all'), prctile(op, 97.5, 'all')])

% set(gca,'ColorScale','log')

% hold on; plot(variedParam(1).range, exp(variedParam(1).range/1.1)-1); plot(variedParam(1).range, exp((variedParam(1).range-1)*5)+15);


%% Plot with better colormap

logPvalData = log10( squeeze(op(1,:,:))');

figure; 
imagesc(xParamvec, yParamvec, logPvalData, 'AlphaData', ~isnan(squeeze(op(1,:,:))'))
set(gca,'YDir','normal')
cb = colorbar(); cb.Label.String = cbLabel2;
xlabel(xName,'Interpreter','none')
ylabel(yName,'Interpreter','none')
title(analysisTitle)

N = 256; n = N/2;
cm = NaN(N,3);
cm(:,1) = [ones(n,1);linspace(1,0,N-n)';];
cm(:,2) = [linspace(0,1,n)';linspace(1,0,N-n)']; 
cm(:,3) = [linspace(0,1,n)';ones(N-n,1)]; 


alpha = 0.05;
% bonferroniCorr = false
% if bonferroniCorr; alpha = 0.05./numel(logPvalData); end

set(gca,'clim',[log10(alpha)*2 0])
set(gcf,'colormap',cm)
cb = colorbar('Direction','reverse','Ticks',[log10(alpha/10),log10(alpha),log10(alpha*10)],'TickLabels',[alpha/10,alpha,alpha*10]);
cb.Label.String = 'KS test p-value';
title(''); xlabel('Mean cluster membership'); ylabel('Number of clusters')

% set(gca, 'XTick',xParamvec, 'XTickLabel',num2str(xParamvec')) % For grid starting at mnc=1.5

%% Plot slice of above p-val matrix

% ithYval = 1:numel(yParamvec); 
%  ithYval = [2, 4,  8,  12]; 
ithYval = [2]; 
xx = xParamvec;
yy = squeeze(op(1,:,ithYval))';
figure; hold on; 
plot(xx, yy, '-o'); yline(0.05)
xlabel('Mean cluster membership'); ylabel('KS test p-value');
set(gca, 'YScale', 'log')
legend({num2str( yParamvec(ithYval)')}, 'Location', 'Best')

%% Extra parameter grid plots

figure; 
imagesc(xParamvec, yParamvec, mean(allNetDecodeSlopes, 3)', 'AlphaData', ~isnan(mean(allNetDecodeSlopes, 3)'))
set(gca,'YDir','normal')
cb = colorbar(); cb.Label.String = 'Mean abs. decode slope';
xlabel(xName,'Interpreter','none'); ylabel(yName,'Interpreter','none')


figure; 
imagesc(xParamvec, yParamvec, mean(allNetDecodeSlopes_zscored, 3)', 'AlphaData', ~isnan(mean(allNetDecodeSlopes_zscored, 3)'))
set(gca,'YDir','normal')
cb = colorbar(); cb.Label.String = 'Shuffle scored decode slope';
xlabel(xName,'Interpreter','none'); ylabel(yName,'Interpreter','none')


figure; 
imagesc(xParamvec, yParamvec, mean(allNetDecodeEntropy, 3)', 'AlphaData', ~isnan(mean(allNetDecodeEntropy, 3)'))
set(gca,'YDir','normal')
cb = colorbar(); cb.Label.String = 'Mean entropy';
xlabel(xName,'Interpreter','none'); ylabel(yName,'Interpreter','none')


figure; 
imagesc(xParamvec, yParamvec, mean(allNetDecodeVar, 3)', 'AlphaData', ~isnan(mean(allNetDecodeVar, 3)'))
set(gca,'YDir','normal')
cb = colorbar(); cb.Label.String = 'Mean decode variance';
xlabel(xName,'Interpreter','none'); ylabel(yName,'Interpreter','none')


%% Plot net-wise scatter, if single parameter point was run above


% figure; scatterhist(allNetMedianDiff(:), allNetPvals(:), 'Group', allNetGroupID(:), 'Kernel','on')
figure; scatter(allNetMedianDiff(:), allNetPvals(:), [], allNetGroupID(:))
title('All networks from all parameter points')
xlabel('median(actual R^2)-median(shuff R^2)')
ylabel('KStest p-value')

%{
figure; scatterhist(allNetKSstat(:), allNetPvals(:), 'Group', allNetGroupID(:), 'Kernel','on')
xlabel('KS-statistic')
ylabel('KStest p-value')
%}

%% Plot net-wise scatter, if multiple parameter points were run

try
    [a, b] = find(squeeze(op(1,:,:))'<0.05); ind1 = [a]; ind2 = [b];
    ind1 = [2]; ind2 = [4];

    X = squeeze(allNetMedianDiff(ind1,ind2,:));
    Y = squeeze(allNetPvals(ind1,ind2,:)); 
    ID = squeeze(allNetGroupID(ind1,ind2,:));
    figure; scatterhist(X(:), Y(:), 'Group', ID(:), 'Kernel','on')
    % figure; scatter(X(:), Y(:), [], ID(:))
    title(['Nets from parameter point index (', num2str(ind1), ', ', num2str(ind2), ')'])
    xlabel('median(actual R^2)-median(shuff R^2)')
    ylabel('KStest p-value')
end


