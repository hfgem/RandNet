

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


minPeakRate = 2; % minimum peak PF rate to be considered a place cell


xParamvec = variedParam(1).range;
xName = variedParam(1).name;
yParamvec = variedParam(2).range;
yName = variedParam(2).name;


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
        parfor ithNet = 1:size(resultsStruct, 3)
            
            % Get matrix of PFs
            try
                PFmat = PFresultsStruct(ithParam1, ithParam2, ithNet).results{1}.linfields;
                E_indices = PFresultsStruct(ithParam1, ithParam2, ithNet).results{1}.E_indices;
            catch
                disp(['PF data error ', num2str(ithParam1), num2str(ithParam2), num2str(ithNet)])
                PFmat = [];
            end
            

            % Calculate place field score
            if ~isempty(PFmat)
                
                % Reformat matrix of PFs to struct
                network = struct; 
                linfields = {};
                network.E_indices = E_indices;
                network.all_indices = 1:parameters.n;
                day = 1; epoch = 1; tetrode = 1; tr = 1;
                linfields{day}{epoch}{tetrode}{1}{tr}(:,1) = pfsim.gridxvals*100; % convert to cm
                for ithCell = network.E_indices
                    linfields{day}{epoch}{tetrode}{ithCell}{tr}(:,5) = PFmat(ithCell,:);
                end
            
                PFscore = calculate_linfieldsScore(linfields, pfsim, pfsim, network);
                
                fprintf('PF score %0.4f %1.0f %1.0f %1.0f \n', PFscore, ithParam1, ithParam2, ithNet)
                
                temp(1, ithNet) = PFscore;
                %temp(2, ithNet) = nan;

            else
                temp(1, ithNet) = nan;
                %temp(2, ithNet) = nan;
            end
            
        end
        
        %{
        temp(1,:)
        nanmean(temp(1,:))
        min(temp(1,:))
        keyboard
        %}
        
        op(1, ithParam1, ithParam2) = nanmean(temp(1,:));
        op(2, ithParam1, ithParam2) = min(temp(1,:));
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
