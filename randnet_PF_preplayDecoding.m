%% Script to run Bayesian decoding on preplay sim from randnet_PF
% This script can be run after randnet_PF to perform Bayesian decoding on
% the preplay sequences

saveFlag = 1; % 1 to save simulation results
selectPath = 0; % 1 to select save destination, 0 to save in current dir
plotResults = 1; % 1 to plot basic simulation results

if saveFlag & selectPath
    save_path = uigetdir('/Users/hannahgermaine/Documents/PhD/','Select Save Folder'); %Have user input where they'd like the output stored
else
    save_path = [pwd, '/results/randnet_PF/DecodingExamples/'];
end

addpath(genpath('functions'))


%% Set up data structures needed for Jadhav lab code, from randnet_PF.m sim

day = 1; eprun = 2; tet = 1;

tetinfo{day}{eprun}{tet}.numcells = parameters.n_E;

if numel(linfields{day})<eprun
    linfields{day}{eprun} = linfields{day}{1};
end
if numel(linfields{day}{eprun}{:})==parameters.n
    linfields{day}{eprun}{tet}(network.I_indices) = [];
end

ripple{day}{eprun}.starttime = trialResults.events(:,1).*parameters.dt;
ripple{day}{eprun}.endtime = trialResults.events(:,2).*parameters.dt;

for cellid = 1:numel(network.E_indices)
    % spikes{day}{eprun}{tet}{cellid}.data(:,1) = V_m(network.E_indices(cellid),:)>= parameters.V_th; % spike times
    spikes{day}{eprun}{tet}{cellid}.data(:,1) = find(V_m(network.E_indices(cellid),:)>= parameters.V_th).*parameters.dt; % spike times
end


%% Save above data structures with appropriate names

animalprefix = datestr(startTime,'yyyy-mm-ddTHH-MM')
if ~exist([save_path, animalprefix, '_direct/'])
    mkdir([save_path, animalprefix, '_direct/'])
end

%save([savedir, '/cellinfo'], cellinfo)
% 01 indicates day one
save([save_path, animalprefix, '_direct/', animalprefix, 'tetinfo.mat'], 'tetinfo')
save([save_path, animalprefix, '_direct/', animalprefix, 'linfields01.mat'], 'linfields')
save([save_path, animalprefix, '_direct/', animalprefix, 'rippletime01.mat'], 'ripple')
save([save_path, animalprefix, '_direct/', animalprefix, 'spikes01.mat'], 'spikes')
% save([save_path, animalprefix, '_direct/', animalprefix, 'params.mat'], 'parameters', 'pfsim')


%% Run preplay decoding analysis
% animalprefix = '2022-06-29T12-40'

day=1; ep=2;

cellcountthresh = 5; % at least 5 cells fired as a candidate event
savedata = saveFlag % save data = 1; not save = 0
figopt = 1 % 1 = plot decoding result for each event; 0 = no figure generated
shuffleIterations = 500 % 1500 for standard decoding

 savedata = 0 % save data = 1; not save = 0

warning off
replaytrajectory = preplay_decoding_CA1_singleday(animalprefix,day,ep,cellcountthresh, save_path, savedata, figopt, shuffleIterations);
warning on


%% Plot preplay decoding results

%{
animalprefix = '2022-06-28T14-09';
save_path = [pwd, '/results/randnet_PF/DecodingExamples/', animalprefix, '_direct/'];% directory for saving results
load([save_path, animalprefix, 'replaydecode_CA1_01_02'])
%}

day=1; ep=2;

pvals = replaytrajectory{day}{ep}.pvalue(:,1);
figure; histogram(pvals, 10)

frame = 1; track = 1;
% replaytrajectory{day}{ep}.shuffle_rsquare{frame}{track}

allshuff_rvals = vertcat(replaytrajectory{day}{ep}.shuffle_rsquare{:});
allshuff_rvals = allshuff_rvals(:,1); % take just forward traj

rvals_preplay = replaytrajectory{day}{ep}.rsquare(:,1);
rvals_shuffle = vertcat(allshuff_rvals{:,1});
figure; hold on; 
ecdf(rvals_preplay)
ecdf(rvals_shuffle)

[H,P,KSSTAT] = kstest2(rvals_preplay, rvals_shuffle)
legend({'Preplays', 'Shuffles'}, 'Location', 'Best')


%% rval and max jump threshold pval matrix

allshuff_rvals = vertcat(replaytrajectory{day}{ep}.shuffle_rsquare{:});
allshuff_rvals = allshuff_rvals(:,1); % take just forward traj

allshuff_jumps = vertcat(replaytrajectory{day}{ep}.shuffle_maxJump{:});
allshuff_jumps = allshuff_jumps(:,1); % take just forward traj

rvals_preplay = replaytrajectory{day}{ep}.rsquare(:,1);
jump_preplay = replaytrajectory{day}{ep}.maxJump(:,1);

rvalThresh_vec = 0:0.1:1;
jumpThres_vec = 0:0.1:1;

op = zeros(numel(jumpThres_vec), numel(rvalThresh_vec));
for ij = 1:numel(jumpThres_vec)
    
    for ir = 1:numel(rvalThresh_vec)
        nActPass = mean( [rvals_preplay>rvalThresh_vec(ir)] & [jump_preplay<jumpThres_vec(ij)]);
        nShuffPass = zeros(1, size(allshuff_jumps{1}, 2));
        for ithShuf = 1:size(allshuff_jumps{1}, 2)
            nShuff_jump = cellfun(@(x) [x(ithShuf)<jumpThres_vec(ij)], allshuff_jumps);
            nShuff_rval = cellfun(@(x) [x(ithShuf)>rvalThresh_vec(ir)], allshuff_rvals);
            nShuffPass(ithShuf) = mean( nShuff_jump & nShuff_rval );
        end

        op(ij, ir) = 1 - mean(nActPass>nShuffPass);
        if [sum(nShuffPass)==0] & [sum(nActPass)==0]
            op(ij, ir) = nan;
        end
    end
    ij
end


figure; imagesc(rvalThresh_vec, jumpThres_vec, op', 'AlphaData', ~isnan(op')); colorbar
xlabel('max jump')
ylabel('r^2')

% Plot with better colormap
figure; 
imagesc(rvalThresh_vec, jumpThres_vec, log10(op'), 'AlphaData', ~isnan(op'))
% set(gca,'YDir','normal')
cb = colorbar(); %cb.Label.String = cbLabel2;
xlabel('max jump')
ylabel('r^2')
%title(analysisTitle)

N = 256; n = N/2;
cm = NaN(N,3);
cm(:,1) = [ones(n,1);linspace(1,0,N-n)';];
cm(:,2) = [linspace(0,1,n)';linspace(1,0,N-n)']; 
cm(:,3) = [linspace(0,1,n)';ones(N-n,1)]; 

set(gca,'clim',[log10(.05)*2 0])
set(gcf,'colormap',cm)
colorbar
colorbar('Direction','reverse','Ticks',[log10(.005),log10(.05),log10(.5)],'TickLabels',[.005,.05,.5])

hold on; 
[xnan, ynan] = find(isnan(op));
scatter(rvalThresh_vec(xnan), jumpThres_vec(ynan), 300, 'x')

