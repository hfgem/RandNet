

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


xParamvec = variedParam(1).range;
xName = variedParam(1).name;
yParamvec = variedParam(2).range;
yName = variedParam(2).name;

analysisTitle = 'PF score';
cbLabel1 = 'Mean PF score';
cbLabel2 = 'Best PF score';

op = nan(2, numel(xParamvec), numel(yParamvec));

figure(515); hold on
imagesc(xParamvec, yParamvec, squeeze(op(1,:,:))', 'AlphaData', ~isnan(squeeze(op(1,:,:))')); drawnow
set(gca,'YDir','normal')
cb = colorbar();
xlabel(xName,'Interpreter','none')
ylabel(yName,'Interpreter','none')

rng(1)
gcp
tic
for ithParam1 = 1:size(resultsStruct, 1)
    for ithParam2 = 1:size(resultsStruct, 2)
        
        temp = nan(1, num_nets);
        for ithNet = 1:size(resultsStruct, 3)
            
            if ~isempty(resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.pMat)
                
                keyboard
                %pvals = resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.pvalue(:,1);
                %figure; histogram(pvals, 10)

                allshuff_rvals = vertcat(resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.shuffle_rsquare{:});
                allshuff_rvals = allshuff_rvals(:,1); % take just forward traj

                rvals_preplay = resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.rsquare(:,1);
                rvals_shuffle = vertcat(allshuff_rvals{:,1});
                figure; hold on; ecdf(rvals_preplay); ecdf(rvals_shuffle); legend({'Preplays', 'Shuffles'}, 'Location', 'Best')

                [H,P,KSSTAT] = kstest2(rvals_preplay, rvals_shuffle);

                temp(1, ithNet) = P;
                temp(2, ithNet) = KSSTAT;

            else
                temp(1, ithNet) = nan;
                temp(2, ithNet) = KSSTAT;
            end
                
        end

        op(1, ithParam1, ithParam2) = nanmean(temp(1,:));
        op(2, ithParam1, ithParam2) = nanmean(temp(2,:));
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


