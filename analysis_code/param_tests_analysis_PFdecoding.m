
if isfolder('functions')
    addpath('functions')
else
    addpath( fullfile( '..', 'functions' ) )
end


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

weightedDecodes = 1;

%downSample = 0
%downSampleFracEvents = 0.5;

xParamvec = variedParam(1).range;
xName = variedParam(1).name;
yParamvec = variedParam(2).range;
yName = variedParam(2).name;

if weightedDecodes
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
            if [isfield(resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory, 'pMat') && ...
                ~isempty(resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.pMat)]

                % keyboard
                %pvals = resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.pvalue(:,1);
                %figure; histogram(pvals, 10)

                if weightedDecodes
                    allshuff_rvals = vertcat(resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.shuffle_weightedR2{:});
                    allshuff_rvals = allshuff_rvals(:,1); % take just forward traj
                    rvals_shuffle = vertcat([allshuff_rvals{:,1}]');
                    rvals_preplay = resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.weightedR2(:,1);
                else
                    allshuff_rvals = vertcat(resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.shuffle_rsquare{:});
                    allshuff_rvals = allshuff_rvals(:,1); % take just forward traj
                    rvals_shuffle = vertcat(allshuff_rvals{:,1});
                    rvals_preplay = resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.rsquare(:,1);
                end
                
                if plotAllCDFs
                    figure; hold on; ecdf(rvals_preplay); ecdf(rvals_shuffle); legend({'Preplays', 'Shuffles'}, 'Location', 'Best')
                end
                
                [H,P,KSSTAT] = kstest2(rvals_preplay, rvals_shuffle);
                temp(1, ithNet) = P;
                temp(2, ithNet) = KSSTAT;
                
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
                
                pvals_preplay = resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.pvalue(:,1);
                nEvents_accum = nEvents_accum + numel(pvals_preplay);
                nSigEvents_accum = nSigEvents_accum + sum(pvals_preplay<0.05);
                
                %%
                
                slopes_preplay = resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.slopes(:,1);
                % slopes_shuffle = vertcat(resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.shuffle_slopes{:});
                % slopes_shuffle = slopes_shuffle(:,1); % take just forward traj
                allNetDecodeSlopes(ithParam1, ithParam2, ithNet) = mean(abs(slopes_preplay));
                % allNetDecodeSlopes_zscored = 0;

                entropy_preplay = resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.Entropy(:,1);
                allNetDecodeEntropy(ithParam1, ithParam2, ithNet) = mean(entropy_preplay);
                
                temp = 0;
                for ithEvent = 1:numel(resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.pMat)
                    temp = temp + mean(var(resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.pMat{ithEvent}{1}.pMat, [], 1));
                end
                allNetDecodeVar(ithParam1, ithParam2, ithNet) = temp./ numel(resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.pMat);
                                
                
                %%
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


% Plot with better colormap
figure; 
imagesc(xParamvec, yParamvec, log10(squeeze(op(1,:,:))'), 'AlphaData', ~isnan(squeeze(op(1,:,:))'))
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

set(gca,'clim',[log10(.05)*2 0])
set(gcf,'colormap',cm)
colorbar
colorbar('Direction','reverse','Ticks',[log10(.005),log10(.05),log10(.5)],'TickLabels',[.005,.05,.5])


%% Extra parameter grid plots

figure; 
imagesc(xParamvec, yParamvec, mean(allNetDecodeSlopes, 3)', 'AlphaData', ~isnan(mean(allNetDecodeSlopes, 3)'))
set(gca,'YDir','normal')
cb = colorbar(); cb.Label.String = 'Mean abs. decode slope';
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

figure; scatterhist(allNetMedianDiff(:), allNetPvals(:), 'Group', allNetGroupID(:), 'Kernel','on')
xlabel('median(actual R^2)-median(shuff R^2)')
ylabel('KStest p-value')

figure; scatterhist(allNetKSstat(:), allNetPvals(:), 'Group', allNetGroupID(:), 'Kernel','on')
xlabel('KS-statistic')
ylabel('KStest p-value')

%% Plot net-wise scatter, if multiple parameter points were run

try
    [a, b] = find(squeeze(op(1,:,:))'<0.05); ind1 = [a]; ind2 = [b];
    ind1 = [2]; ind2 = [4];

    X = squeeze(allNetMedianDiff(ind1,ind2,:));
    Y = squeeze(allNetPvals(ind1,ind2,:)); 
    ID = squeeze(allNetGroupID(ind1,ind2,:));
    figure; scatterhist(X(:), Y(:), 'Group', ID(:), 'Kernel','on')
    xlabel('median(actual R^2)-median(shuff R^2)')
    ylabel('KStest p-value')
end


