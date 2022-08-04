function replaytrajectory = preplay_decoding_CA1_singleday(animalprefix,day,ep,cellcountthresh, savedir, savedata, figopt, shuffleIterations)
% Adapted from code for Shin et al., 2019, Jadhav Lab
%
% INPUTS:
%    animalprefix = animal prefix.
%    day = experimental day.
%    ep = epoch.
%    cellcountthresh = mininum number of cells active for considering as a
%                      cadidate event, usually = 5
%    savedir = directory for saving results
%    savedata = save results, 1 = save, 0 = not save
%    figopt = plot figures, 1 = plot, 0 = no plot

% This function reguires the files: 
% cellinfo, linfields01, rippletime01, spikes01, and tetinfo
%
% Producing rippletime01 requires the script replay_SWRtime_singleday.m,
% which requires the files: 
% pos01, ripples01, tetinfo
%
% Assumptions:
% assumes all simulation data is on tetrode 1
% Assumes all simulation data does not exclude any cells


%% Analysis parameters

tBinSz = 10; %ms, default temporal bin in ms [typically around 15ms], hard coded
minEventDur = 50; % ms, exclude events shorter than this
wellcutoff = 0; %cm, remove reward-well regions (15cm around start and end); or 0cm without exclusion
minPeakRate = 3; %Hz, minimum peak rate to include cell

downSampleCells = 0;
downSampledFrac = 0.1;

dispFlag = 1; % display event analysis progression

%% set animal directory
dir = [savedir, animalprefix, '_direct/'];
exclude_list = [];

%% 
if ~mod(ep,2)
   eprun = ep;% using current epoch to construct place-field templates for RUN epoch
elseif ep == 1
   eprun = 2;% using the first run epoch for pre-task sleep epoch
else
   eprun = ep-1;% using the previous run epoch for sleep epoch (post-task sleep)
end

%%
%-----match neurons across epochs-----%

switch animalprefix
    case {'ER1','KL8','JS14','JS15','JS17','JS21'} % if non-sim, then find cells across ep
        [~, hpidx] = matchidx_acrossep_singleday(dir, animalprefix, day,exclude_list); %(tet, cell)
    otherwise % if sim data, all cells are correct
        %load(sprintf('%s%stetinfo.mat',dir,animalprefix)); % get num tets, get num cells
        tetinfo = load(sprintf('%s%stetinfo.mat',dir,animalprefix), 'tetinfo'); % get num tets, get num cells
        tetinfo = tetinfo.tetinfo; 
        numTets = numel(tetinfo{1}{2});
        
        if ~downSampleCells % keep all cells
            numCells = tetinfo{1}{2}{1}.numcells;
            hpidx = [repmat(numTets, numCells, 1), [1:numCells]'];
        else % Use only fraction of cells for decoding
            numCells = round(downSampledFrac * tetinfo{1}{2}{1}.numcells);
            hpidx = [repmat(numTets, numCells, 1), [1:numCells]'];
        end
end

hpnum = length(hpidx(:,1));

%%
%-----create the ratemaps [nPosBin x nHPCells]-----%
rm = []; % ratemap matrix
pm = []; % position matrix
tm = []; % track matrix
nTracks = 1;
cellidxm = [];
%load(sprintf('%s%slinfields0%d.mat',dir,animalprefix,day)); % get linearized place fields
linfields = load(sprintf('%s%slinfields0%d.mat',dir,animalprefix,day), 'linfields'); % get num tets, get num cells
linfields = linfields.linfields; 
for i = 1:hpnum
      cind = hpidx(i,:);
      if (length(linfields{day}{eprun})>= cind(1))
            if  (length(linfields{day}{eprun}{cind(1)})>= cind(2))
                linfield1 = linfields{day}{eprun}{cind(1)}{cind(2)};
            else 
                linfield1 =[];
            end
      else
            linfield1=[];
      end
      
      if ~isempty(linfield1)
           linfield_hp = [];
           lintrack_hp = [];
           pos_hp = [];
           for track = 1:nTracks
                temp1 = linfield1{track};
                pos1 = temp1(:,1);
                lintrack1 = ones(size(pos1))*track;
                occnormrate1 = temp1(:,5);
                linfield_hp = [linfield_hp;occnormrate1];
                pos_hp = [pos_hp;pos1];
                lintrack_hp = [lintrack_hp;lintrack1];
           end
           if (max(linfield_hp) >= minPeakRate) % peak firing rate max larger than 3 Hz
               a = find(isnan(linfield_hp));
               %pad nan
               if ~isempty(a)
                    [lo,hi]= findcontiguous(a);  %find contiguous NaNs
                    for ii = 1:length(lo) 
                        if lo(ii) > 1 & hi(ii) < length(linfield_hp)
                            fill = linspace(linfield_hp(lo(ii)-1), ...
                                linfield_hp(hi(ii)+1), hi(ii)-lo(ii)+1);
                            linfield_hp(lo(ii):hi(ii)) = fill;
                        end
                    end
               end
               rm = [rm;linfield_hp'];
               pm = [pm;pos_hp'];
               tm = [tm;lintrack_hp'];
               cellidxm = [cellidxm; cind];
           end
      end
end
rm = rm'; %[nPosBin x nHPCells]
pm = pm';
tm = tm';

% remove reward-well regions (15cm around start and end)
for i = 1:nTracks
    pm_traj = pm(find(tm == i));
    maxpos = max(max(pm_traj));
    rm(find(tm == i & pm <= wellcutoff)) = 0;
    rm(find(tm == i & pm >= maxpos-wellcutoff)) = 0;
end

rm = rm+ (eps.^8); %Add a small number so there are no zeros
expecSpk =rm.*tBinSz./1000; %[nPos x nCells] Expected number of spikes per bin
hpnum = length(rm(1,:)); % update cell number (cell with < 3Hz peak rate excluded)


%%
%-----create the event matrix during SWRs-----%
spikes = loaddatastruct(dir, animalprefix, 'spikes', day); % get spikes


% get ripple time
% load(sprintf('%s%srippletime0%d.mat',dir,animalprefix,day));
ripple = load(sprintf('%s%srippletime0%d.mat',dir,animalprefix,day), 'ripple'); % get num tets, get num cells
ripple = ripple.ripple; 

rip = ripple{day}{ep}; 
if numel(rip.starttime)>numel(rip.endtime)
    riptimes(:,1) = rip.starttime(1:numel(rip.endtime));
    riptimes(:,2) = rip.endtime;    
else
    riptimes(:,1) = rip.starttime;
    riptimes(:,2) = rip.endtime;   
end
rip_starttime = 1000*riptimes(:,1);  % in ms

dur = 1000*(riptimes(:,2) - riptimes(:,1)); % event duration
keepidx = find(dur >= minEventDur);%at least 5 bins, 50 ms for 10ms bins; exclude events < 50ms
rip_starttime = rip_starttime(keepidx);
riptimes = riptimes(keepidx,:);


%%
% loop
if ~isempty(riptimes)
    celldata = [];
    spikecounts = [];
    % cell loop, measure active cells during each event
    for cellcount = 1:hpnum
        index = [day,ep,cellidxm(cellcount,:)] ;
        if ~isempty(spikes{index(1)}{index(2)}{index(3)}{index(4)}.data)
            spiketimes = spikes{index(1)}{index(2)}{index(3)}{index(4)}.data(:,1);
        else
            spiketimes = [];
        end
        spikebins = periodAssign(spiketimes, riptimes(:,[1 2]));
        if ~isempty(spiketimes)
            validspikes = find(spikebins);
            spiketimes = spiketimes(validspikes);
            spikebins = spikebins(validspikes);
            tmpcelldata = [spiketimes spikebins];
        end
        if ~isempty(spiketimes)
            tmpcelldata(:,3) = cellcount;
        else 
            tmpcelldata = [0 0 cellcount];
        end
        celldata = [celldata; tmpcelldata];
        spikecount = zeros(1,size(riptimes,1));
        for i = 1:length(spikebins)
            spikecount(spikebins(i)) = spikecount(spikebins(i))+1;
        end
        spikecounts = [spikecounts; spikecount];
    end
    cellcounts = sum((spikecounts > 0));
    eventindex = find(cellcounts >= cellcountthresh); % exclude events < 5cells active
    
    % event loop
    for event = 1:length(eventindex)
        
        % show current event number
        if dispFlag
            disp(['Event ', num2str(event), ' of ', num2str(length(eventindex))])
        end
        
        cellsi = celldata(find(celldata(:,2)==eventindex(event)),3); 
        [cellsi,ia] = unique(cellsi,'first');
        [~,sortorder] = sort(ia);
        event_cellSeq = cellsi(sortorder);
        tmpind = find(celldata(:,2) == eventindex(event));
        spiketimes = celldata(tmpind,1);
        cellindex = celldata(tmpind,3);
        %-----create the event matrix during SWRs (spkT{cells}.spiketimes) -----%
        for cell = event_cellSeq'
            validspikeidx = find(cellindex == cell);
            spkT{cell} = spiketimes(validspikeidx).*1000;
        end     
        %--------calculate the posterior probability on an event-by-event basis--------%
        startevent = riptimes(eventindex(event),1).*1000;
        endevent = riptimes(eventindex(event),2).*1000;
        timebins = startevent:tBinSz:endevent; % timebins are the binedges
        nTBin = length(timebins)-1;
        nCell = hpnum;
        spkPerBin = zeros(1,nTBin, nCell); % keep the inactive cells as 0s.
        for nn  = 1:hpnum
            cellInd = nn; %current cell
            if length(spkT) >= cellInd
                if ~isempty(spkT{cellInd})
                    temp = histc(spkT{cellInd}, timebins); %[1 x nTBin x nCell]
                    spkPerBin(1,:,cellInd) = temp(1:end-1);
                end
            end
        end
        nSpkPerTBin = squeeze(sum(spkPerBin,3)); %[nTBin x 1] number of spikes in tBin  
        nPBin = size(rm,1); %N positional bin
        
        % Baye's rules
        expecSpk  = reshape(expecSpk, [nPBin,1, nCell]); %[nPos x 1 x nCell]
        expon = exp(-expecSpk); %Exponent of equation.
        factSpkPerBin = factorial(spkPerBin); %Factorial to divide by

        wrking = bsxfun(@power, expecSpk, spkPerBin); %[nPos x nTbin x nCell]
        wrking = bsxfun(@rdivide, wrking, factSpkPerBin); %[nPos x nTbin x nCell]
        wrking = bsxfun(@times,wrking, expon); %[nPos x nTbin x nCell]
        post = prod(wrking,3); %Non normalised prob [nPos x Tbin]
        
        %{
        figure; imagesc(exp(sum(log(wrking),3))); colorbar
        figure; imagesc(sum(log(wrking),3)); colorbar
        
        logPost = sum( log(wrking) , 3); 
        figure; imagesc(logPost./min(logPost, [], 1)); colorbar
        figure; imagesc(exp( logPost./min(logPost, [], 1) )); colorbar
        logNormPost = exp( logPost./-min(logPost, [], 1) );
        %logNormPost = logNormPost(2:end-1,:);
        figure; imagesc(logNormPost./sum(logNormPost, 1)); colorbar
        tic
        a = prod(vpa(wrking) ,3); toc
        figure; imagesc(double(a./sum(a, 1)))
        
        b =  exp(vpa(sum(log(wrking) ,3)));
        figure; imagesc(double(a./sum(a, 1)))

        keyboard
        %}

        post(:,nSpkPerTBin==0)  =0; % so the posterior matrix can be smoothed.
        post(isnan(post)) = 0;   % remove NaN
                
        trajinfo = mean(tm,2);% trajactory number
        
        if figopt==1
            figure(9342) % plot result
        elseif figopt==2
            myPlotSettings(8, 2)
            figure; sgtitle([animalprefix, ', Day: ', num2str(day), ', Ep: ', num2str(ep), ', Event:', num2str(event)]); 
        end
        % trajectory loop (4 trajectory types in a W-maze)
        for traj = 1:4 
            trajidx = find(trajinfo == traj);
            pMat = post(trajidx,:);% create a posterior matrix for each traj
            distvector = pm(trajidx,1)';

            szPM1 = size(pMat,1);
            szPM2 = size(pMat,2);
            for i = 1:szPM2; if (sum(pMat(:,i))>0)
            pMat(:,i) = pMat(:,i)./sum(pMat(:,i));end;end % normalized across positions to 1 for each time bin 
            nonzerobins = find(nSpkPerTBin > 0);
            
            % Monte Carlo simulation for linear regression
            rvalues = [];
            slopes = [];
            entropy = [];
            interc = [];
            fitpVals = [];
            totalsamples = 10000;
            for rloop = 1:500
                tBinPicks = distsample(totalsamples,nSpkPerTBin);
                regressdata = [];
                for i = 1:length(nonzerobins)
                    if (nSpkPerTBin(nonzerobins(i)) > 0)
                       tmpnumsamples = sum(tBinPicks == nonzerobins(i));
                       if ~isempty(find(pMat(:,nonzerobins(i)) ~= 0))
                           distpicks = distvector(distsample(tmpnumsamples,pMat(:,nonzerobins(i))))';
                           entropy_loop(i) = -nansum((hist(distpicks,0:5:200)./length(distpicks)).*log(hist(distpicks,0:5:200)./length(distpicks)));
                           distpicks(:,2) = nonzerobins(i);
                           regressdata = [regressdata; distpicks];
                       end
                    end
                end
                regressdata(:,3) = 1;

                [b,bint,r,rint,stats] = regress(regressdata(:,1),[regressdata(:,3),regressdata(:,2)]);
                rvalues = [rvalues; stats(1)];
                fitpVals = [fitpVals; stats(3)];
                slopes = [slopes; b(2)];
                interc = [interc; b(1)];
                
                if ~exist('entropy_loop') % JB, added 3/4/21
                   entropy_loop = 0;
                end
                
                entropy = [entropy; mean(entropy_loop)];
            end
                        
            [~, PeakDecodeInd] = max(pMat, [], 1); 
            if ~isempty(PeakDecodeInd)
                maxJump(event,traj) = max(abs(diff(PeakDecodeInd))./size(pMat, 1));
            else
                maxJump(event,traj) = nan;
            end
            
            Res(event,traj) =  mean(rvalues);
            Fit_pVal(event, traj) = mean(fitpVals);
            YInterc(event,traj) = mean(interc);
            Spd(event,traj) = mean(slopes);
            Entropy(event,traj) = mean(entropy);
            pMat_cell{event}{traj}.pMat = pMat;
            pMat_cell{event}{traj}.timevec = 1:szPM2;
            pMat_cell{event}{traj}.posvec = distvector;
            pMat_cell{event}{traj}.timebinsz = tBinSz;
            
            [X,y] = meshgrid([1:size(pMat,2)],[1:size(pMat,1)]); w = pMat;
            mdl = fitlm(X(:),y(:),'Weights',w(:));
            weightedSlope(event,traj) = mdl.Coefficients.Estimate(2);
            weightedR2(event,traj) = mdl.Rsquared.Ordinary;
            
            if figopt
                subplot(1,4,traj)
                imagesc(1:szPM2,distvector,pMat);
                colormap(jet)
                hold on
                plot([1, szPM2], [0, Spd(event,traj) * (szPM2-1)] + (YInterc(event,traj) +1), 'w','linewidth',2);
                hold off
    %             pause(0.1)
                caxis([0 0.1])
            end
            
            %-------Shuffling to get the pvalue for each traj------%
            permbins = nonzerobins;
            srvalues = [];
            smaxJumps = [];
            sslopes = [];
            sweightedSlope = [];
            sweightedR2 = [];
            for iteration = 1:shuffleIterations % 1500
                permbins = permbins(randperm(length(permbins)));% temporal shuffle
                % calculate shuffled pMat
                tmpspkPerBin = zeros(size(spkPerBin));
                tmpspkPerBin(:,permbins,:) = spkPerBin(:,nonzerobins,:);
                tmpfactSpkPerBin = factorial(tmpspkPerBin); %Factorial to divide by
                tmpnSpkPerTBin = squeeze(sum(tmpspkPerBin,3)); %[nTBin x 1] number of spikes in tBin  

                tmpexpecSpk = expecSpk(trajidx,:,:);
                tmpexpon = expon(trajidx,:,:);
                wrking = bsxfun(@power, tmpexpecSpk, tmpspkPerBin); %[nPos x nTbin x nCell]
                wrking = bsxfun(@rdivide, wrking, tmpfactSpkPerBin); %[nPos x nTbin x nCell]
                wrking = bsxfun(@times,wrking, tmpexpon); %[nPos x nTbin x nCell]
                tmppMat = prod(wrking,3); %Non normalised prob [nPos x Tbin]
                tmppMat(:,tmpnSpkPerTBin==0)  =0; % so the posterior matrix can be smoothed.
                tmppMat(isnan(tmppMat)) = 0;  
                for i = 1:szPM2; if (sum(tmppMat(:,i))>0)
                tmppMat(:,i) = tmppMat(:,i)./sum(tmppMat(:,i));end;end % normalized across positions to 1 for each time bin 
                clear wrking tmpfactSpkPerBin tmpexpon
                
                % Monte Carlo simulation for linear regression
                tBinPicks = distsample(totalsamples,nSpkPerTBin);
                regressdata = [];
                for i = 1:length(permbins)
                    if (nSpkPerTBin(nonzerobins(i)) > 0)
                       tmpnumsamples = sum(tBinPicks == nonzerobins(i));
                       if ~isempty(find(pMat(:,nonzerobins(i)) ~= 0))
                           distpicks = distvector(distsample(tmpnumsamples,tmppMat(:,permbins(i))))';
                           distpicks(:,2) = permbins(i);
                           regressdata = [regressdata; distpicks];
                       end
                    end
                end
                regressdata(:,3) = 1;
                [b,bint,r,rint,stats] = regress(regressdata(:,1),[regressdata(:,3),regressdata(:,2)]);
                srvalues = [srvalues; stats(1)];
                sslopes = [sslopes; b(2)];
                
                [~, PeakDecodeInd] = max(tmppMat, [], 1); 
                currentShuffMaxJump = max(abs(diff(PeakDecodeInd))./size(tmppMat, 1));
                smaxJumps = [smaxJumps, currentShuffMaxJump];
                
                [X,y] = meshgrid([1:size(tmppMat,2)],[1:size(tmppMat,1)]); w = tmppMat;
                mdl = fitlm(X(:),y(:),'Weights',w(:));
                sweightedSlope = [sweightedSlope,  mdl.Coefficients.Estimate(2)];
                sweightedR2 = [sweightedR2, mdl.Rsquared.Ordinary];
                
            end
            
            % calculate p-value
            pvalue(event,traj) = sum(Res(event,traj) < srvalues)/length(srvalues);
            
            % Save shuffle stats
            shuffle_Spd{event}{traj} = sslopes;
            shuffle_rvalues{event}{traj} = srvalues;
            shuffle_maxJump{event}{traj} = smaxJumps;
            shuffle_weightedSlope{event}{traj} = sweightedSlope;
            shuffle_weightedR2{event}{traj} = sweightedR2;
        end
        
        [minP,tidx] = min(pvalue(event,:));% find the minimun pvalue
        if minP < 0.05 % significant
            decode_traj(event) = tidx; % trajectory with the minimun pvalue represented
        else
            decode_traj(event) = 0;% no significant traj
        end
        % save cell info during the event
        activecell{event} = cellsi;
        activecellidx{event} = cellidxm(cellsi,:);
        clear wrking factSpkPerBin expon
        
        if figopt==2
            for subplotind = 1:4
                subplot(1,4,subplotind);
                xlabel('Time (10 ms bin)');
                title(num2str(pvalue(event, subplotind)));
                if subplotind == 1
                    ylabel('Position (cm)');
                end
            end
        end
               
    end
end

%%
% structure result
if ~exist('pMat_cell')
    replaytraj.pMat = [];
    replaytraj.maxJump = [];
    replaytraj.rsquare = [];
    replaytraj.FitpVal = [];
    replaytraj.slopes = [];
    replaytraj.YInterc = [];
    replaytraj.Entropy = [];
    replaytraj.eventidx = [];
    replaytraj.besttraj = [];
    replaytraj.pvalue = [];
    replaytraj.activecell = [];
    replaytraj.activecellidx = [];
    replaytraj.sigeventprc = [];
    replaytraj.sigeventnum = [];
    replaytraj.candeventnum = [];
    replaytraj.wellcutoff = wellcutoff;
    replaytraj.tBinSz = tBinSz;
    replaytraj.cellcountthresh = cellcountthresh;
    
    replaytraj.shuffle_rsquare = [];
    replaytraj.shuffle_maxJump = [];
else
    replaytraj.pMat = pMat_cell;
    replaytraj.maxJump = maxJump; % maximum jump in peak decoded pos. between time bins, in frac of track
    replaytraj.rsquare = Res;
    replaytraj.FitpVal = Fit_pVal;
    replaytraj.slopes = Spd;
    replaytraj.YInterc = YInterc;
    replaytraj.Entropy = Entropy;
    replaytraj.eventidx = eventindex;
    replaytraj.besttraj = decode_traj;
    replaytraj.pvalue = pvalue;
    replaytraj.activecell = activecell;
    replaytraj.activecellidx = activecellidx;
    replaytraj.sigeventprc = length(find(decode_traj~=0))./length(decode_traj);
    replaytraj.sigeventnum = length(find(decode_traj~=0));
    replaytraj.candeventnum = length(decode_traj);
    replaytraj.wellcutoff = wellcutoff;
    replaytraj.tBinSz = tBinSz;
    replaytraj.cellcountthresh = cellcountthresh;
    
    replaytraj.weightedSlope = weightedSlope;
    replaytraj.weightedR2 = weightedR2;
    
    replaytraj.shuffle_slopes = shuffle_Spd;
    replaytraj.shuffle_rsquare = shuffle_rvalues;
    replaytraj.shuffle_maxJump = shuffle_maxJump;
    
    replaytraj.shuffle_weightedSlope = shuffle_weightedSlope;
    replaytraj.shuffle_weightedR2 = shuffle_weightedR2;
end

replaytrajectory{day}{ep} = replaytraj;


%%
%---save data ---%
if savedata
   save(sprintf('%s%sreplaydecode_CA1_%02d_%02d.mat', dir,animalprefix,day,ep), 'replaytrajectory');
end