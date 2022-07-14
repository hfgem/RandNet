

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
%downSample = 0
%downSampleFracEvents = 0.5;

xParamvec = variedParam(1).range;
xName = variedParam(1).name;
yParamvec = variedParam(2).range;
yName = variedParam(2).name;

analysisTitle = 'Decoding KS-test';
if combineNetData
    cbLabel1 = 'combined p-val';
    cbLabel2 = 'combined KS-stat';
else
    cbLabel1 = 'Median p-val';
    cbLabel2 = 'Median KS-stat';
end

op = nan(2, numel(xParamvec), numel(yParamvec));

figure(515); hold on
imagesc(xParamvec, yParamvec, squeeze(op(1,:,:))', 'AlphaData', ~isnan(squeeze(op(1,:,:))')); drawnow
set(gca,'YDir','normal')
cb = colorbar();
xlabel(xName,'Interpreter','none')
ylabel(yName,'Interpreter','none')

rng(1)
% gcp
tic
for ithParam1 = 1:size(resultsStruct, 1)
    for ithParam2 = 1:size(resultsStruct, 2)
        
        temp = nan(1, num_nets);
        allRvecs = [];
        allRvecs_shuff = [];
        for ithNet = 1:size(resultsStruct, 3)
            
            % if ~isempty(resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.pMat)
            if [isfield(resultsStruct(ithParam1, ithParam2, ithNet).results, 'replaytrajectory') && ...
                ~isempty(resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.pMat)]

                % keyboard
                %pvals = resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.pvalue(:,1);
                %figure; histogram(pvals, 10)

                allshuff_rvals = vertcat(resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.shuffle_rsquare{:});
                allshuff_rvals = allshuff_rvals(:,1); % take just forward traj

                rvals_preplay = resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.rsquare(:,1);
                rvals_shuffle = vertcat(allshuff_rvals{:,1});
                
                if plotAllCDFs
                    figure; hold on; ecdf(rvals_preplay); ecdf(rvals_shuffle); legend({'Preplays', 'Shuffles'}, 'Location', 'Best')
                end
                
                [H,P,KSSTAT] = kstest2(rvals_preplay, rvals_shuffle);
                temp(1, ithNet) = P;
                temp(2, ithNet) = KSSTAT;
                
                
                allRvecs = [allRvecs, rvals_preplay'];
                allRvecs_shuff = [allRvecs_shuff, rvals_shuffle'];

            else
                temp(1, ithNet) = nan;
                temp(2, ithNet) = KSSTAT;
            end
                
        end
        
        if combineNetData
            if ~isempty(allRvecs)
                allRvecs = allRvecs(1:min(numel(allRvecs), maxNEvents));
                
                [~,p_kstest,KSSTAT] = kstest2(allRvecs, allRvecs_shuff);
                op(1, ithParam1, ithParam2) = p_kstest;
                op(2, ithParam1, ithParam2) = KSSTAT;
            else
                p_kstest = nan;
                op(1, ithParam1, ithParam2) = nan;
                op(2, ithParam1, ithParam2) = nan;
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
