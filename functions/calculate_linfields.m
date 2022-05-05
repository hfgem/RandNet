function [linfields, PFpeaksSequence] = calculate_linfields(opS, parameters, sim, network, plotPFs)
% Calculate linear place fields, based on Jadhav lab methods
%
% Adapted from ReplayNet code, with assumption of simulating a single
% trajectory for a single epoch


%% Set up

% Default indexes, to ensure consistency with data structure used by Jadhav lab
day=1; epoch=1; tetrode=1; 

xPos = sim.xPos; 

winSize = sim.linFieldGaussSD/sim.spatialBin*sim.winScale; % window is 5x the STDev
linTrajBins = numel(sim.gridxvals); % ~40 for arms, ~20 for half of top width, 2 cm bins

%% Calculate occupancy for each trajectory (independent of cell)
accumOcc = zeros(linTrajBins, 1);        % sum of occupancy timesteps along the 4 linear trajectories
for i = 1:numel(sim.t)
    currTraj = 1;
    switch currTraj
        case 1
            [~, linI] = min(( abs( [sim.gridxvals] - [xPos(i)] ) ));
            accumOcc(linI, currTraj) = accumOcc(linI,currTraj) + 1;  
    end

end
accumOcc = accumOcc.*parameters.dt*sim.nTrials;

% figure; plot(accumOcc); legend({'tr1', 'tr2', 'tr3', 'tr4'}); xlabel('Distance from center arm (m)'); ylabel('Un-norm Lin. Place field')
% figure; plot(smoothdata(accumOcc,'gaussian',winSize)); legend({'tr1', 'tr2', 'tr3', 'tr4'}); xlabel('Distance from center arm (m)'); ylabel('Un-norm Lin. Place field')
% keyboard 

%% For each cell, calculate its linear place field
for ithCell = 1:parameters.n %increment through cells

    %spikeTimes = spikes{day}{epoch}{tetrode}{ithCell}.data(:,1);
    spikeTimes = find(squeeze(sum(opS(ithCell,:,:), 3)))*parameters.dt;

    accumTrajSpikes = zeros(linTrajBins, 1); % sum of ithCell's spikes along the 4 linear trajectories

    for i = 1:numel(spikeTimes)
        curX = xPos(round(spikeTimes(i)/parameters.dt)); % xPos at spike time
        [xL, linXind] = min(abs(sim.gridxvals - curX)); % linear xPos index
        currTraj = 1;
        switch currTraj
            case 0
                % stationary
            case 1
                [~, linI] = min(abs( [sim.gridxvals] - [curX] ));
                accumTrajSpikes(linI, currTraj) = accumTrajSpikes(linI,currTraj) +1;
        end
    end
    % figure; scatter(xPos(round(spikeTimes/parameters.dt)), yPos(round(spikeTimes/parameters.dt)));
    % figure; plot(accumTrajSpikes); legend({'tr1', 'tr2', 'tr3', 'tr4'}); xlabel('Distance from center arm (m)'); ylabel('Un-norm Lin. Place field')
    % figure; plot(smoothdata(accumTrajSpikes,'gaussian',winSize)); legend({'tr1', 'tr2', 'tr3', 'tr4'}); xlabel('Distance from center arm (m)'); ylabel('Un-norm Lin. Place field')

    % Calculate linfields values
    for tr = 1
        linfields{day}{epoch}{tetrode}{ithCell}{tr}(:,1) = [0:sim.spatialBin:(linTrajBins*sim.spatialBin-sim.spatialBin)]' *100; % (cm) linear distance from center well 
        linfields{day}{epoch}{tetrode}{ithCell}{tr}(:,2) = accumOcc(:,tr); % occupancy
        linfields{day}{epoch}{tetrode}{ithCell}{tr}(:,3) = accumTrajSpikes(:,tr); % spike count
        linfields{day}{epoch}{tetrode}{ithCell}{tr}(:,4) = accumTrajSpikes(:,tr)./accumOcc(:,tr); % occupancy-normalized firing rate
        linfields{day}{epoch}{tetrode}{ithCell}{tr}(:,5) = smoothdata(accumTrajSpikes(:,tr),'gaussian',winSize)./smoothdata(accumOcc(:,tr),'gaussian',winSize); % smoothed occupancy-normalized firing rate
        linfields{day}{epoch}{tetrode}{ithCell}{tr}(:,6) = smoothdata(accumOcc(:,tr),'gaussian',winSize); % smoothed occupancy
        linfields{day}{epoch}{tetrode}{ithCell}{tr}(:,7) = smoothdata(accumTrajSpikes(:,tr),'gaussian',winSize); % smoothed spike count
    end

    %{
    tr = 1;
    figure; plot(linfields{day}{epoch}{tetrode}{ithCell}{tr}(:,1), linfields{day}{epoch}{tetrode}{ithCell}{tr}(:,2))
    figure; plot(linfields{day}{epoch}{tetrode}{ithCell}{tr}(:,1), linfields{day}{epoch}{tetrode}{ithCell}{tr}(:,3))
    figure; plot(linfields{day}{epoch}{tetrode}{ithCell}{tr}(:,1), linfields{day}{epoch}{tetrode}{ithCell}{tr}(:,4))
    figure; plot(linfields{day}{epoch}{tetrode}{ithCell}{tr}(:,1), linfields{day}{epoch}{tetrode}{ithCell}{tr}(:,5))
    figure; plot(linfields{day}{epoch}{tetrode}{ithCell}{tr}(:,1), linfields{day}{epoch}{tetrode}{ithCell}{tr}(:,6))
    figure; plot(linfields{day}{epoch}{tetrode}{ithCell}{tr}(:,1), linfields{day}{epoch}{tetrode}{ithCell}{tr}(:,7))
    keyboard
    %}

end


% Extract place field order from linfields struct
opTemp = [];
for ithCell = 1:parameters.n
    PF = linfields{day}{epoch}{tetrode}{ithCell}{tr}(:,5);
    %if sum(PF>0)
        opTemp= [opTemp;PF'];
    %end
end
PFmat_E = opTemp(network.E_indices,:);

row_all_zeros1 = find(all( PFmat_E==0, 2)) ;
row_n_all_zeros1 = find(~all( PFmat_E==0, 2)) ;

[peakRate,peakRateLocation] = max(squeeze(PFmat_E(row_n_all_zeros1,:)), [], 2);
[B,sortedCellIndsbyPeakRateLocation] = sort(peakRateLocation, 'descend');
PFpeaksSequence = [row_n_all_zeros1(sortedCellIndsbyPeakRateLocation); row_all_zeros1];



normRates = 1;
if normRates
    rateDenom1 = max(PFmat_E(PFpeaksSequence,:), [], 2);
    caxmax = 1;
else
    rateDenom1 = 1;
    caxmax = max(PFmat_E, [], 'all');
end
if plotPFs
    figure; imagesc( PFmat_E(PFpeaksSequence,:)./rateDenom1 ); title('Env1, sort Env1'); colorbar; caxis([0, caxmax])
    xlabel('Position (2 cm bin)'); ylabel('Cell (sorted))');
    
    figure; histogram(rateDenom1); 
    xlabel('Peak PF rate (Hz)'); ylabel('E cells (count)');
end

peakRate = max(PFmat_E(PFpeaksSequence,:), [], 2);
PFpeaksSequence(peakRate<sim.minPeakRate) = nan;

end
