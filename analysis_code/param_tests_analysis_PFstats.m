

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

op = nan(4, numel(xParamvec), numel(yParamvec));

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
        
        temp = nan(4, num_nets);
        for ithNet = 1:size(resultsStruct, 3)
            
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
                
                % Select only active, excitatory cells
                PFinds = [ ismember(1:size(PFmat, 1), E_indices)' & (max(PFmat, [], 2)>=minPeakRate) ];
                ePFmat = PFmat(PFinds,:);

                
                % 'MeanRate'
                meanPeakRate= mean( max(ePFmat, [], 2) );
                fprintf('PF mean Peak Rate %0.4f %1.0f %1.0f %1.0f \n', meanPeakRate, ithParam1, ithParam2, ithNet)
                temp(1, ithNet) = meanPeakRate;  
                    
                % 'kstest'
                accum_ksstat = 0;
                for i = 1:size(ePFmat, 1)
                    [~,p,ksstat,~] = kstest( ( ePFmat(i,:)-mean(ePFmat(i,:), 2) )./std(ePFmat(i,:), [], 2) );
                    accum_ksstat = accum_ksstat + ksstat;
                end
                %fprintf('PF KS-stat %0.4f %1.0f %1.0f %1.0f \n', accum_ksstat/size(ePFmat, 1), ithParam1, ithParam2, ithNet)
                temp(2, ithNet) = accum_ksstat/size(ePFmat, 1);                
                    
                % 'sparsity'                
                cellSparsity = mean( ePFmat<=[0.25*max(ePFmat, [], 2)], 2 );
                %fprintf('PF sparsity %0.4f %1.0f %1.0f %1.0f \n', mean(cellSparsity), ithParam1, ithParam2, ithNet)
                temp(3, ithNet) = mean(cellSparsity);

                % 'information'
                spatialInfo = nanmean( [ePFmat./mean(ePFmat, 2)] .* log(( ePFmat+eps )./mean(ePFmat, 2) ), 2 );
                %fprintf('PF info. %0.4f %1.0f %1.0f %1.0f \n', nanmean(spatialInfo), ithParam1, ithParam2, ithNet)
                temp(4, ithNet) = nanmean(spatialInfo);                

            else
                temp(:, ithNet) = nans(4,1);
            end
            
        end
        
        %{
        temp(1,:)
        nanmean(temp(1,:))
        min(temp(1,:))
        keyboard
        %}
        
        op(:, ithParam1, ithParam2) = nanmean(temp, 2);

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
cb = colorbar(); cb.Label.String = '';
xlabel(xName,'Interpreter','none')
ylabel(yName,'Interpreter','none')

figure; 
imagesc(xParamvec, yParamvec, squeeze(op(2,:,:))', 'AlphaData', ~isnan(squeeze(op(2,:,:))'))
set(gca,'YDir','normal')
cb = colorbar(); cb.Label.String = '';
xlabel(xName,'Interpreter','none')
ylabel(yName,'Interpreter','none')

figure; 
imagesc(xParamvec, yParamvec, squeeze(op(3,:,:))', 'AlphaData', ~isnan(squeeze(op(2,:,:))'))
set(gca,'YDir','normal')
cb = colorbar(); cb.Label.String = '';
xlabel(xName,'Interpreter','none')
ylabel(yName,'Interpreter','none')

figure; 
imagesc(xParamvec, yParamvec, squeeze(op(4,:,:))', 'AlphaData', ~isnan(squeeze(op(2,:,:))'))
set(gca,'YDir','normal')
cb = colorbar(); cb.Label.String = '';
xlabel(xName,'Interpreter','none')
ylabel(yName,'Interpreter','none')


%{

figure;
subplot(2,2,1)
imagesc(variedParam(paramPlot1).range, variedParam(paramPlot2).range, frac_partic', 'AlphaData', ~isnan(frac_partic'))
set(gca,'YDir','normal')
c1 = colorbar(); c1.Label.String = 'Fraction of neurons';
title('Frac. firing (event)'); xlabel(variedParam(paramPlot1).name,'Interpreter','none'); ylabel(variedParam(paramPlot2).name,'Interpreter','none')

subplot(2,2,2)
imagesc(variedParam(paramPlot1).range, variedParam(paramPlot2).range, avg_fr', 'AlphaData', ~isnan(avg_fr'))
set(gca,'YDir','normal')
c2 = colorbar(); c2.Label.String = "Hz";
title('Mean spike rate (trial)'); xlabel(variedParam(paramPlot1).name,'Interpreter','none'); ylabel(variedParam(paramPlot2).name,'Interpreter','none')

subplot(2,2,3)
imagesc(variedParam(paramPlot1).range, variedParam(paramPlot2).range, avg_event_length', 'AlphaData', ~isnan(avg_event_length'))
set(gca,'YDir','normal')
c3 = colorbar(); c3.Label.String = "Seconds";
title('Mean Event Length'); xlabel(variedParam(paramPlot1).name,'Interpreter','none'); ylabel(variedParam(paramPlot2).name,'Interpreter','none')
clear c1 c2 c3

subplot(2,2,4)
imagesc(variedParam(paramPlot1).range, variedParam(paramPlot2).range, avg_n_events'/parameters.t_max, 'AlphaData', ~isnan(avg_n_events'))
set(gca,'YDir','normal')
c3 = colorbar(); c3.Label.String = "nEvents / s";
title('Mean event frequency'); xlabel(variedParam(paramPlot1).name,'Interpreter','none'); ylabel(variedParam(paramPlot2).name,'Interpreter','none')
clear c1 c2 c3

%}

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
