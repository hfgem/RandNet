function linfields = calculate_linfields(opS, parameters, sim)
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

end
