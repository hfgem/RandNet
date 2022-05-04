function PFScore = calculate_linfieldsScore(linfields, parameters, sim, network)

%% For each cells' place field, calculate KS-statistic (KSop) and other statistics
day = 1; epoch = 1; tetrode = 1; tr = 1;

% Unconstrained gaussian leads to high, persistant rates
% Constrain gaussian fit parameters to realistic place field properties
% 3 < peak firing rate < 30 Hz
% 0 < mean firing rate position < 100 cm
% 0 < standard deviation < 10 ~= firing above 25% of max covers ~25% of track

% gaus = @(x, a1, b1, c1,vo) a1*exp(-((x-b1)/c1).^2);
% figure; plot(posBins, gaus(posBins, 100, 50, 10) ); hold on; yline(25)
% a1 = amplitude = peak firing rate 3<peak<30
% b1 = mean = location of peak 0<mean<100 cm
% c1 = standard deviation = width of place field STD<10

%sim.mainSeed
%sim.gaussFO.Method

gaussFO = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',sim.gaussFOLower,...
               'Upper',sim.gaussFOUpper...
               ); % [peak amplitude, position of peak on track, standard deviation of peak]

           
posBins = linfields{day}{epoch}{tetrode}{1}{tr}(:,1);
    
KSop = []; % Quantify similarity to gaussian
maxFR = [];
sparsity = [];
minFR = [];
spatialInfo = [];
allPFs = [];
gaussRsqrs = [];
for ithCell = 1:parameters.n
    PF = linfields{day}{epoch}{tetrode}{ithCell}{tr}(:,5);
    
    if std(PF)~=0
        [h,p,ksstat,cv] = kstest((PF-mean(PF))/std(PF));
    else % cell is silent
        ksstat = 0;
    end
    
    KSop(ithCell) = ksstat;
    maxFR(ithCell) = max(PF);
    minFR(ithCell) = min(PF);
    sparsity(ithCell) = mean(PF<=(0.25*max(PF)));
    
    if mean(PF)==0
        spatialInfo(ithCell) = 0;
    else
        spatialInfo(ithCell) = mean( [PF/mean(PF)].*log((PF+eps)/mean(PF)) );
        %figure; plot(PF); title(num2str(spatialInfo(ithCell)))
    end
    allPFs(:,ithCell)= PF;
    
    [fitobject,gof,output] = fit(posBins, PF, 'gauss1', gaussFO); % up to gauss8, for sum of 8 gaussians
    
    gaussRsqrs(ithCell) = gof.rsquare;
        
end

% figure; imagesc(allPFs)
% figure; scatter(KSop, sparsity)
% figure; histogram(KSop); title('Distribition of KS-statistics')
% figure; scatter(KSop, dipOP); xlabel('KS-statitic'); ylabel('dip statistic')



%% Calculate objective score, based on PF similarities to Gaussian

% Set all improper gaussRsqrs to minimum
gaussRsqrs(isnan(gaussRsqrs))= -1;
gaussRsqrs(isinf(gaussRsqrs))= -1;
gaussRsqrs(gaussRsqrs<-1) = -1;

% If a cell never fires, make it nan to exclude from analysis
maxFR(maxFR==0)=nan;

PFScore = 0;
for i = 1:numel(posBins)
    eCellInds = ismember(network.all_indices, network.E_indices);
    Inds = logical( ismember(network.all_indices, network.E_indices).* [allPFs(i,:)==maxFR] ); % get index of cells whose peak is at spatial bin i 

    locMax = max(gaussRsqrs(logical(Inds))); % Determine best gaussian correlation of out of Inds cells

    if ~isempty(locMax)
        MeasureOfhighFiring = (1- [abs( maxFR(Inds)-sim.peakTarget ) / sim.peakTarget]); %index for proximity of a peak of 15 Hz
        PFScore = PFScore + -max( gaussRsqrs(Inds)+MeasureOfhighFiring ); % accumulate to the score the max gaussRsqr out of cells with peak at position i
    end

end
PFScore = PFScore/numel(posBins);
% keyboard

end