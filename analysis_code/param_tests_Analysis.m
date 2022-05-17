
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

parameters.n
parameters.del_G_sra
variedParam(1).name
variedParam(2).name

correlationType = 'Pearson'; 
nShuffMultiplier = 3;
minNShuff = 50; % minimum nShuffles
minDetectedSequences = 5

useMaxSeq = 0
maxDetectedSequences = 100

xParamvec = variedParam(1).range;
xName = variedParam(1).name;
yParamvec = variedParam(2).range;
yName = variedParam(2).name;

tic
op = nan(2, numel(xParamvec), numel(yParamvec));

figure(515); hold on
imagesc(xParamvec, yParamvec, squeeze(op(1,:,:))', 'AlphaData', ~isnan(squeeze(op(1,:,:))')); drawnow
set(gca,'YDir','normal')
cb = colorbar();
xlabel(xName,'Interpreter','none')
ylabel(yName,'Interpreter','none')

rng(1)
for ithParam1 = 1:size(resultsStruct, 1)
    
    for ithParam2 = 1:size(resultsStruct, 2)
        
        temp = nan(2, num_nets);
        for ithNet = 1:size(resultsStruct, 3)
            % keyboard
            
            try
            mat = resultsStruct(ithParam1, ithParam2, ithNet).results.ranksVec;
            x = mat./max(mat);
            catch
                % resultsStruct(ithParam1, ithParam2, ithNet).results
                x = [];
            end
            
            if size(x, 2)>=minDetectedSequences

                
                % Sequence-by-sequence correlation analysis
                cbLabel1 = 'p-val';
                cbLabel2 = 'KS-stat';
                analysisTitle = 'KS-test, against shuffle';
                
                correlationType = 'Pearson'; % Pearson, Kendall, Spearman
                if useMaxSeq
                    if size(x, 2) > maxDetectedSequences
                        x = x(:, randsample(size(x, 2), maxDetectedSequences) );
                    end
                end
                nSeq = size(x, 2); % number of sequences in x
                rMat_full = corr(x, 'rows','pairwise');
                rMat = tril(rMat_full, -1); rMat(rMat==0)=nan;
                %{
                rMat = nan(nSeq, nSeq);
                for i = 1:(nSeq-1)
                    for j = i+1:nSeq
                        [r,p] = corr(x(:,i),x(:,j),'type',correlationType, 'rows','complete');
                        rMat(i, j) = r;
                        rMat(j, i) = r;
                    end
                end
                %}
                % figure; imagesc(rMat); colorbar
                
                % Shuffled sequence-by-sequence correlation analysi
                nShuf = max(nShuffMultiplier * size(x, 2), minNShuff);
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
                
                rMat_full_shuff = corr(x_shuff, 'rows','pairwise');
                rMat_shuff = tril(rMat_full_shuff, -1); rMat_shuff(rMat_shuff==0)=nan;
                %{
                rMat_shuff = nan(nShuf, nShuf);
                for i = 1:(nShuf-1)
                    for j = i+1:nShuf
                        [r,p] = corr(x_shuff(:,i),x_shuff(:,j),'type',correlationType, 'rows','complete');
                        % shuffledRmat(i, j) = r;
                        rMat_shuff(j, i) = r;
                    end
                end
                %}
                % figure; imagesc(rMat_shuff); colorbar
                
                % figure; hold on; ecdf(rMat_shuff(:),'Bounds','on'); cdfplot(rMat(:)); 
                % legend({'Shuffle', 'Actual'}, 'Location', 'Best')
                
                if sum(~isnan(rMat), 'all')>0
                    [~,p_kstest,KSSTAT] = kstest2(rMat(:), rMat_shuff(:));
                else
                    p_kstest = nan;
                    KSSTAT = nan;
                end
                temp(1, ithNet) = p_kstest;
                temp(2, ithNet) = KSSTAT;

                %{               
                % Mean Rel Rank linear correlation
                analysisTitle = 'Mean rel. rank corr.';
                [vals, inds] = sort(nanmean(x, 2), 'ascend');
                indsRep = repmat(inds, 1, size(x, 2));
                xSort = x(inds,:);
                xVals = x(inds,:);
                yVals = repmat(1:size(x,1), 1, size(x,2));
                % figure; scatter(xVals(:), yVals(:)); hold on; scatter(nanmean(xVals, 2), 1:size(xVals, 1)); xlabel('Rel. rank'); ylabel('Neuron')
                
                mdl = fitlm(xVals(:), yVals(:)); % fieldnames(mdl)
                temp(ithNet) = mdl.Rsquared.adjusted; cbLabel='r^2 adj.';% mdl.Rsquared.adjusted, mdl.LogLikelihood, mdl.Coefficients.pValue(2)
                %temp(ithNet) = mdl.LogLikelihood; cbLabel='LogLikelihood';% mdl.Rsquared.adjusted, mdl.LogLikelihood, mdl.Coefficients.pValue(2)
                %temp(ithNet) = mdl.Coefficients.pValue(2); cbLabel='p-val';% mdl.Rsquared.adjusted, mdl.LogLikelihood, mdl.Coefficients.pValue(2)
                %}
                
                %{
                % Mean Rel Rank linear correlation
                analysisTitle = 'Shuffle: Mean rel. rank corr.';
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
                [vals, inds] = sort(nanmean(x_shuff, 2), 'ascend');
                indsRep = repmat(inds, 1, size(x_shuff, 2));
                xSort = x_shuff(inds,:);
                xVals = x_shuff(inds,:);
                yVals = repmat(1:size(x_shuff,1), 1, size(x_shuff,2));
                % figure; scatter(xVals(:), yVals(:)); hold on; scatter(nanmean(xVals, 2), 1:size(xVals, 1)); xlabel('Rel. rank'); ylabel('Neuron')
                
                mdl = fitlm(xVals(:), yVals(:)); % fieldnames(mdl)
                %temp(ithNet) = mdl.Rsquared.adjusted; cbLabel='r^2 adj.';% mdl.Rsquared.adjusted, mdl.LogLikelihood, mdl.Coefficients.pValue(2)
                temp(ithNet) = mdl.LogLikelihood; cbLabel='LogLikelihood';% mdl.Rsquared.adjusted, mdl.LogLikelihood, mdl.Coefficients.pValue(2)
                %temp(ithNet) = mdl.Coefficients.pValue(2); cbLabel='p-val';% mdl.Rsquared.adjusted, mdl.LogLikelihood, mdl.Coefficients.pValue(2)
                %}
                                
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
                temp(1, ithNet) = nan;
                temp(2, ithNet) = nan;
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


%%
X = squeeze(op(1,:,:))';
figure; plot(X); set(gca,'YScale','log')

%% Detect threshold crossing, then plot fitted exponential curve

thresh = 0.4;
x = squeeze(op(2,:,:))>thresh; 
% figure; imagesc(mncVec, nClustersVec, x'); set(gca, 'YDir', 'normal')
% figure; imagesc(mncVec, nClustersVec,  diff(x, [], 2)'); set(gca, 'YDir', 'normal')
[xvalinds, yvalinds] = find( diff(x, [], 2) )

ft = fittype('a*exp(b*x) + c');
% [fexp_offset, fexp_offset_G] = fit(mncVec(xvalinds)', nClustersVec(yvalinds)', ft);
[fexp, fexp_G] = fit(xParamvec(xvalinds)', yParamvec(yvalinds)','exp1');

hold on; scatter(xParamvec(xvalinds), yParamvec(yvalinds), 'r')
plot(xParamvec, fexp(xParamvec), 'r')

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