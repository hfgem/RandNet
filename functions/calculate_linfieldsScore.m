function PFScore = calculate_linfieldsScore(linfields, parameters, sim, network, varargin)

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

%% Default parameters:
plotPFs = 0;

%% Read in optional parameters, to overwrite above defaults
for i=1:2:length(varargin)
    switch varargin{i}
        case 'plotPFs'
            plotPFs = varargin{i+1};
        otherwise
            error('plot_ixj_sequences: Unknown input')
    end
end


%% Fit each PF to a constrained Gaussian

% Defaults for TolFun and TolX are 1e-6
gaussFO = fitoptions('Method','NonlinearLeastSquares', ...
                    'Lower',sim.gaussFOLower, ...
                    'Upper',sim.gaussFOUpper, ...
                    'TolFun', 1e-3, ...
                    'TolX', 1e-3 ...
                    ); % [peak amplitude, position of peak on track, standard deviation of peak]

           
posBins = linfields{day}{epoch}{tetrode}{1}{tr}(:,1);
    
maxFR = nan(1, parameters.n);
allPFs = nan(numel(posBins), parameters.n);
gaussRsqrs = nan(1, parameters.n);
for ithCell = network.E_indices
    PF = linfields{day}{epoch}{tetrode}{ithCell}{tr}(:,5);
    allPFs(:,ithCell)= PF;
    maxFR(ithCell) = max(PF);
    
	[~,gof,~] = fit(posBins, PF, 'gauss1', gaussFO); % up to gauss8, for sum of 8 gaussians
	gaussRsqrs(ithCell) = gof.rsquare;
        
end


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

%% Plot place field
if plotPFs
        
    PFmat = allPFs';
    PFmat_E = PFmat(network.E_indices,:);

    row_all_zeros1 = find(all( PFmat_E==0, 2)) ;
    row_n_all_zeros1 = find(~all( PFmat_E==0, 2)) ;

    [peakRate,peakRateLocation] = max(squeeze(PFmat_E(row_n_all_zeros1,:)), [], 2);
    [B,sortedCellIndsbyPeakRateLocation] = sort(peakRateLocation, 'descend');
    PFpeaksSequence = [row_n_all_zeros1(sortedCellIndsbyPeakRateLocation); row_all_zeros1];

    [peakRate, peakRateLocation_all] = max(PFmat_E, [], 2);


    normRates = 1;
    if normRates
        rateDenom1 = max(PFmat_E(PFpeaksSequence,:), [], 2);
        caxmax = 1;
    else
        rateDenom1 = 1;
        caxmax = max(PFmat_E, [], 'all');
    end

    figure; histogram(rateDenom1); 
    xlabel('Peak PF rate (Hz)'); ylabel('E cells (count)');

    figure; imagesc( PFmat_E(PFpeaksSequence,:)./rateDenom1 ); 
    title(['Score=', num2str(PFScore)]); 
    colorbar; caxis([0, caxmax])
    xlabel('Position (2 cm bin)'); ylabel('Cell (sorted)');
end

end