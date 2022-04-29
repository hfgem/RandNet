
% Temp code,
% addpath('C:\Users\Jordan\Box\Data\Replay project\RandNet\ParameterSweeps')
% load('results_2022-04-13T19-46.mat') % mnc X nClusters
% load('results_2022-04-15T01-42.mat') % Gin X Wee

% load('results_2022-04-28T19-37.mat')
% load('results_2022-04-28T23-08.mat')

%% Analyze sequence properties from resultsStruct
%
% resultsStruct(param1Ind, param2Ind, ithNet)
% 
% num_nets 
% variedParam

correlationType = 'Pearson'; 

xParamvec = variedParam(1).range;
xName = variedParam(1).name;
yParamvec = variedParam(2).range;
yName = variedParam(2).name;

tic
op = zeros(numel(xParamvec), numel(yParamvec));
for ithParam1 = 1:size(resultsStruct, 1)
    
    for ithParam2 = 1:size(resultsStruct, 2)
        
        temp = zeros(num_nets, 1);
        for ithNet = 1:size(resultsStruct, 3)
            % keyboard
            
            try
            mat = resultsStruct(ithParam1, ithParam2, ithNet).results.ranksVec;
            x = mat./max(mat);
            catch
                % resultsStruct(ithParam1, ithParam2, ithNet).results
                x = [];
            end
            
            if size(x, 2)>2

                %keyboard
                analysisTitle = 'Mean rel. rank corr.';
                [vals, inds] = sort(nanmean(x, 2), 'ascend');
                indsRep = repmat(inds, 1, size(x, 2));
                xSort = x(inds,:);
                xVals = x(inds,:);
                yVals = repmat(1:size(x,1), 1, size(x,2));
                % figure; scatter(xVals(:), yVals(:))
                mdl = fitlm(xVals(:), yVals(:)); % fieldnames(mdl)
                temp(ithNet) = mdl.Rsquared.adjusted; % mdl.Rsquared.adjusted, mdl.LogLikelihood, mdl.Coefficients.pValue(2)
                 %[rho,pval] = corr(vals(~isnan(vals)),inds(~isnan(vals)));
                 %temp(ithNet) = pval; % mdl.Rsquared.adjusted, mdl.LogLikelihood, mdl.Coefficients.pValue(2)

                %{
                % Sequence-by-sequence correlation analysis
                correlationType = 'Pearson'; % Pearson, Kendall, Spearman
                nSeq = size(x, 2); % number of sequences in x
                rMat = zeros(nSeq, nSeq);
                for i = 1:nSeq
                    for j = 1:nSeq
                        [r,p] = corr(x(:,i),x(:,j),'type',correlationType, 'rows','complete');
                        rMat(i, j) = r;
                    end
                end
                temp_rMat = rMat;
                temp_rMat(ismembertol(temp_rMat, 1, 10^-12))=nan;
                overallMeanCorr = nanmean(temp_rMat, 'all');
                temp(ithNet) = overallMeanCorr;
                %}
                
                % Shuffled sequence-by-sequence correlation analysis
                %{
                nShuf = 1* size(x, 2);
                x_shuff = zeros( size(x, 1), nShuf);
                for i = 1:nShuf
                    % randomly select from an actual sequence
                    randSeq = x(:, randi(size(x, 2))); 
                    firedInd = find(~isnan(randSeq));
                    % Randomly permute only those cells that actually fired
                    shufSeq = nan(size(randSeq));
                    shufSeq(firedInd) = randSeq(firedInd(randperm(numel(firedInd))));
                    x_shuff(:,i) = shufSeq;
                end
                shuffledRmat = zeros(nShuf, nShuf);
                for i = 1:nShuf
                    for j = 1:nShuf
                        [r,p] = corr(x_shuff(:,i),x_shuff(:,j),'type',correlationType, 'rows','complete');
                        shuffledRmat(i, j) = r;
                    end
                end
                temp_rMat = shuffledRmat;
                temp_rMat(ismembertol(temp_rMat, 1, 10^-12))=nan;
                overallMeanCorr = nanmean(temp_rMat, 'all');
                temp(ithNet) = overallMeanCorr;
                %}
                
                %{
                % Dim. Red. clustering analysis
                Y_tsne = tsne(rMat);  
                Y_tsne_norm = (Y_tsne- mean(Y_tsne)) ./ std(Y_tsne) ;
                % figure; scatter(Y_tsne(:,1),Y_tsne(:,2))
                
                if ~isempty(Y_tsne)
                     [~,p,KSSTAT] = kstest(Y_tsne_norm);
                    % [BF, BC] = bimodalitycoeff(Y_tsne_norm);
                    % [BF, BC] = bimodalitycoeff(Y_tsne);
                else
                    p = nan; KSSTAT = nan;
                    BF = nan; BC = nan;
                end
                temp(ithNet) = nanmean(BC);
                %}

            else
                temp(ithNet) = nan;
            end
            
        end
        
        op(ithParam1, ithParam2) = nanmean(temp);
        
        
    end
end
runTime = toc;
disp([ 'Runtime: ', datestr(datenum(0,0,0,0,0,runTime),'HH:MM:SS') ])
% disp( duration(0, 0, runTime) )

figure; 
imagesc(xParamvec, yParamvec, op', 'AlphaData', ~isnan(op'))
set(gca,'YDir','normal')
colorbar
xlabel(xName)
ylabel(yName)
title(analysisTitle)

% caxis([prctile(op, 2.5, 'all'), prctile(op, 97.5, 'all')])

% set(gca,'ColorScale','log')



%% % %%%%%%%%%%%%
% %% Functions %%
 %%% %%%%%%%%%%%%

function [BF, BC] = bimodalitycoeff(x)
% check if x is vector and if it is - 
% represent it as a column-vector
if isvector(x), x = x(:); end

% determine the data size along its first dimension
N = size(x, 1);

% determine the data skewness (unbiased)
S = skewness(x, 0);
% determine the data kurtosis (unbiased)
% (the normal distribution has kurtosis of zero)
K = kurtosis(x, 0) - 3;

% calculate the bimodality coefficient (unbiased)
BC = (S.^2 + 1) ./ (K + 3* ([(N-1)^2] / [(N-2)*(N-3)]) );

% determine the bimodality flag (using +5% margin)
BF = BC > 1.05*5/9;
end