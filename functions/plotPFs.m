function plotPFs(allLinfields, allPFpeaksSeq, network, pfsim)

day = 1; epoch = 1; tetrode = 1; tr = 1;


%% Convert place fields to matrix

allEnvPFs = [];

for ithEnv = 1:pfsim.nEnvironments

    linfields = allLinfields{ithEnv};
    posBins = linfields{day}{epoch}{tetrode}{1}{tr}(:,1);


    maxFR = nan(1, pfsim.n);
    allPFs = nan(numel(posBins), pfsim.n);
    for ithCell = network.E_indices
        PF = linfields{day}{epoch}{tetrode}{ithCell}{tr}(:,5);
        allPFs(:,ithCell)= PF;
        maxFR(ithCell) = max(PF); 
    end

    allEnvPFs{ithEnv} = allPFs;
end

%% Plot place fields by self-peak order

for ithEnv = 1:pfsim.nEnvironments 

PFmat = allEnvPFs{ithEnv}';
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

%figure; histogram(rateDenom1); 
%xlabel('Peak PF rate (Hz)'); ylabel('E cells (count)');

figure; imagesc( PFmat_E(PFpeaksSequence,:)./rateDenom1 ); 
title(['Env ID ', num2str(ithEnv)]); 
colorbar; caxis([0, caxmax])
xlabel('Position (2 cm bin)'); ylabel('Cell (sorted)');


end

%% Plot place fields by env1 peaks

for ithEnv = 1:pfsim.nEnvironments 

PFmat = allEnvPFs{ithEnv}';
PFmat_E = PFmat(network.E_indices,:);

%{
row_all_zeros1 = find(all( PFmat_E==0, 2)) ;
row_n_all_zeros1 = find(~all( PFmat_E==0, 2)) ;

[peakRate,peakRateLocation] = max(squeeze(PFmat_E(row_n_all_zeros1,:)), [], 2);
peakRateLocation = allPFpeaksSeq(:,ithEnv);

[B,sortedCellIndsbyPeakRateLocation] = sort(peakRateLocation, 'descend');
PFpeaksSequence = [row_n_all_zeros1(sortedCellIndsbyPeakRateLocation); row_all_zeros1];
%}

PFpeaksSequence = allPFpeaksSeq(:,1);
[peakRate, peakRateLocation_all] = max(PFmat_E, [], 2);


normRates = 1;
if normRates
    rateDenom1 = max(PFmat_E(PFpeaksSequence,:), [], 2);
    caxmax = 1;
else
    rateDenom1 = 1;
    caxmax = max(PFmat_E, [], 'all');
end

%figure; histogram(rateDenom1); 
%xlabel('Peak PF rate (Hz)'); ylabel('E cells (count)');

figure; imagesc( PFmat_E(PFpeaksSequence,:)./rateDenom1 ); 
title(['Env ID ', num2str(ithEnv)]); 
colorbar; caxis([0, caxmax])
xlabel('Position (2 cm bin)'); ylabel('Cell (sorted)');


end

%% Plot place fields by all different peak sequences

figure
for ithEnvSequence = 1:pfsim.nEnvironments 
    
    PFpeaksSequence = allPFpeaksSeq(:,ithEnvSequence);
    if mod(ithEnvSequence,2)==0
        PFpeaksSequence = flip(PFpeaksSequence);
    end

    for ithEnvData = 1:pfsim.nEnvironments 

        PFmat = allEnvPFs{ithEnvData}';
        PFmat_E = PFmat(network.E_indices,:);

        [peakRate, peakRateLocation_all] = max(PFmat_E, [], 2);


        normRates = 1;
        if normRates
            rateDenom1 = max(PFmat_E(PFpeaksSequence,:), [], 2);
            caxmax = 1;
        else
            rateDenom1 = 1;
            caxmax = max(PFmat_E, [], 'all');
        end

        %figure; histogram(rateDenom1); 
        %xlabel('Peak PF rate (Hz)'); ylabel('E cells (count)');

        ithPlot = sub2ind( [pfsim.nEnvironments, pfsim.nEnvironments], ithEnvSequence, ithEnvData);
        subplot(pfsim.nEnvironments, pfsim.nEnvironments, ithPlot);
        imagesc( PFmat_E(PFpeaksSequence,:)./rateDenom1 ); 
        % colorbar; caxis([0, caxmax])
        %title(['Env ID ', num2str(ithEnv)]); 
        %xlabel('Position (2 cm bin)'); ylabel('Cell (sorted)');


    end
end
    
end