

%% Analyze sequence and place fields properties of results from randnet_PF_param_tests.m
%
% resultsStruct(param1Ind, param2Ind, ithNet)
% PFresultsStruct(param1Ind, param2Ind, ithNet)
%
% num_nets 
% variedParam
%

parameters.n
parameters.del_G_sra
variedParam(1).name
variedParam(2).name

correlationType = 'Pearson'; 
nShuffMultiplier = 3; % 
minNShuff = 50; % minimum nShuffles, if nSequences*nShuffMultiplier is lower than minNShuff
minDetectedSequences = 5; % minimum preplay sequences to run analyis

minPeakRate = 2; % minimum peak PF rate to be considered a place cell

useMeanPFDensity = 0
useMaxSeq = 0
maxDetectedSequences = 100;

shuffleMethod = 3 % 1 to shuffle sequences, 2 to shuffle PFs, 3 to compare each preplay sequence to many shuffled PFs
nPFshuffles = 100;
sigAlpha = 0.05;

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
            
            try
                mat = resultsStruct(ithParam1, ithParam2, ithNet).results.ranksVec;
                x = mat./max(mat);
            catch
                x = [];
            end
            
            % Get PF sequence, either based on peak or expected location
            try
                PFmat = PFresultsStruct(ithParam1, ithParam2, ithNet).results{1}.linfields;
                E_indices = PFresultsStruct(ithParam1, ithParam2, ithNet).results{1}.E_indices;
                
                [M,I] = max(PFmat(E_indices,:), [], 2); 

                if useMeanPFDensity % Sort by location of mean PF density
                    
                    PFexpectedLocation = sum( PFmat./sum(PFmat, 2) .*([1:size(PFmat, 2)]), 2) ;
                    PFexpectedLocation = PFexpectedLocation(E_indices);
                    PFexpectedLocation(M<minPeakRate)=nan
                    % figure; histogram(PFexpectedLocation)
                    PFseq = PFexpectedLocation;
                    
                else % Sort by location of PF peak
                    I(M<minPeakRate) = nan;
                    [Isorted, PFseq] = sort(I); 
                end
                
                % figure; imagesc(PFmat(E_indices(PFseq),:)./M(PFseq));
                
                PFseq_Rank = PFseq;
                PFseq_Rank(isnan(Isorted))=nan; 
                PFseq_Rank = PFseq_Rank./max(PFseq_Rank);

            catch
                PFseq = [];
            end
            
            % Calculate correlation to preplay sequences
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
                rVec = corr(PFseq_Rank, x, 'rows', 'pairwise');

                % Shuffle either preplay sequences of PF sequence
                if shuffleMethod==1 % Compare PF sequence to shuffled preplays
                    % Shuffled sequence-by-sequence correlation analysis
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
                    rVec_shuff = corr(PFseq_Rank, x_shuff, 'rows', 'pairwise');


                    if sum(~isnan(rVec), 'all')>0
                        % [~,p_kstest,KSSTAT] = kstest2(rMat(:), rMat_shuff(:));
                        [~,p_kstest,KSSTAT] = kstest2(rVec, rVec_shuff);
                    else % Compare preplay sequences to shuffled PF sequences
                        p_kstest = nan;
                        KSSTAT = nan;
                    end
                    
                elseif shuffleMethod==2 % Compare preplay sequences to shuffled PF sequences
                    
                    shuffOP = zeros(1, nPFshuffles);
                    for i = 1:nPFshuffles
                        PF_shuff = PFseq_Rank(randperm(length(PFseq_Rank))); 
                        rVec_shuff = corr(PF_shuff, x, 'rows', 'pairwise');
                        
                        %[~,p_kstest,KSSTAT] = kstest2(rVec, rVec_shuff);
                        %shuffOP(i) = p_kstest;
                        shuffOP(i) = mean(rVec_shuff);
                        
                    end
                    % figure; histogram(shuffOP, 50); xline(mean(rVec))
                    
                    p_kstest =  1 - mean( mean(rVec) > shuffOP );
                    KSSTAT = nan;
                    
                elseif shuffleMethod==3 % Compare individual preplay sequences to shuffled PF sequences
                    
                    for i = 1:nPFshuffles
                        PF_shuffMat(:,i) = PFseq_Rank(randperm(length(PFseq_Rank))); 
                    end
                    
                    sequences_rVec_shuff = zeros(nPFshuffles, size(x, 2));
                    for i = 1:size(x, 2)
                        sequences_rVec_shuff(:,i) = corr(x(:,i), PF_shuffMat, 'rows', 'pairwise');
                        % figure; histogram(sequences_rVec_shuff(:,i)); xline(rVec(i))
                    end
                    sequencePval = 1 - mean(rVec > sequences_rVec_shuff);
                    p_kstest = mean( sequencePval < sigAlpha ) ;
                    KSSTAT = nan;
                    
                end
                
                temp(1, ithNet) = p_kstest;
                temp(2, ithNet) = KSSTAT;
                

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


%% Calculate PF score across 




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