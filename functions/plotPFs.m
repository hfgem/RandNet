function plotPFs(allLinfields, allPFpeaksSeq, network, pfsim)
%
%
% Note: the correlation of the sequence of peaks may be giving falsely high
% correlations because (since there are only 50 spatial PF bins) there were
% many cells that had peaks in the same spatial bin, and tied cells were
% ordered by their index, for all environments
% Should tied cell's be randomized? (by environment? Or by bootstrapping?)


day = 1; epoch = 1; tetrode = 1; tr = 1;

minPeakRate = 0;

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
%{
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
%}


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
        
        mapCorr = nanmean(diag(corr( allEnvPFs{ithEnvSequence},  allEnvPFs{ithEnvData})))
        title(['Map corr=', num2str(mapCorr, '%0.2f')])
        
        if ithEnvData==ithEnvSequence
            nthEnv = ceil(ithEnvData/2);
            if mod(ithEnvData,2)==0; trajDir='leftward'; 
            else; trajDir='rightward'; end
            title(['Env ', num2str(nthEnv), ', ', trajDir])
            
        end
        % colorbar; caxis([0, caxmax])
        %title(['Env ID ', num2str(ithEnv)]); 
        %xlabel('Position (2 cm bin)'); ylabel('Cell (sorted)');
        keyboard
        if ithEnvSequence==1 && ithEnvData==pfsim.nEnvironments
            xlabel('Position (2 cm bin)'); ylabel('Cell (sorted by column''s order)');
        end

    end
end

%% Plot peak sequence correlations
    
figure
for ithEnv1 = 1:pfsim.nEnvironments 
    
    PFpeaksSequence1 = allPFpeaksSeq(:,ithEnv1);
    if mod(ithEnv1,2)==0
        PFpeaksSequence1 = flip(PFpeaksSequence1);
    end
    PFmat = allEnvPFs{ithEnv1}';
    PFmat_E = PFmat(network.E_indices,:);
    [peakRate, peakRateLocation_all] = max(PFmat_E, [], 2);  
    
    PFpeaksSequence1(peakRate<minPeakRate) = nan;
        
    for ithEnv2 = 1:pfsim.nEnvironments 

        PFpeaksSequence2 = allPFpeaksSeq(:,ithEnv2);
        if mod(ithEnv2,2)==0
            PFpeaksSequence2 = flip(PFpeaksSequence2);
        end
        PFmat = allEnvPFs{ithEnv2}';
        PFmat_E = PFmat(network.E_indices,:);
        [peakRate, peakRateLocation_all] = max(PFmat_E, [], 2);  

        PFpeaksSequence2(peakRate<minPeakRate) = nan;
    

        ithPlot = sub2ind( [pfsim.nEnvironments, pfsim.nEnvironments], ithEnv1, ithEnv2);
        subplot(pfsim.nEnvironments, pfsim.nEnvironments, ithPlot);
        scatter(PFpeaksSequence1, PFpeaksSequence2 ); 
        
        [RHO,PVAL] = corr(PFpeaksSequence1, PFpeaksSequence2, 'rows','complete');
        % title([num2str(ithEnv1), num2str(ithEnv2), 'p=', num2str(PVAL)])
        title(['r=', num2str(RHO, '%0.2f')])
        % colorbar; caxis([0, caxmax])
        %title(['Env ID ', num2str(ithEnv)]); 
        %xlabel('Position (2 cm bin)'); ylabel('Cell (sorted)');

        %if ithEnv2==1; ylabel(['Env ', num2str(ithEnv2), ' data']); end

    end
end

end