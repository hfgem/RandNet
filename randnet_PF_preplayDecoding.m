
%% Script to run Bayesian decoding on preplay sim from randnet_PF

% Set up data structures needed for Jadhav lab code, from randnet_PF.m sim

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

% ep = 1;
for cellid = 1:numel(network.E_indices)
    % spikes{day}{eprun}{tet}{cellid}.data(:,1) = V_m(network.E_indices(cellid),:)>= parameters.V_th; % spike times
    spikes{day}{eprun}{tet}{cellid}.data(:,1) = find(V_m(network.E_indices(cellid),:)>= parameters.V_th).*parameters.dt; % spike times
end

%%  Save above data structures with appropriate names
userDir = char(java.lang.System.getProperty('user.home'));
savedir = [userDir, '/Desktop/temp/'];% directory for saving results

animalprefix = 'simy';

%save([savedir, '/cellinfo'], cellinfo)
% 01 indicates day one
save([savedir, animalprefix, '_direct/', animalprefix, 'tetinfo.mat'], 'tetinfo')
save([savedir, animalprefix, '_direct/', animalprefix, 'linfields01.mat'], 'linfields')
save([savedir, animalprefix, '_direct/', animalprefix, 'rippletime01.mat'], 'ripple')
save([savedir, animalprefix, '_direct/', animalprefix, 'spikes01.mat'], 'spikes')


%%
addpath("C:\Users\Jordan\Desktop")

% replay_decoding_CA1_singleday.m


userDir = char(java.lang.System.getProperty('user.home'));
savedir = [userDir, '/Desktop/temp/simy_direct/x'];% directory for saving results


animalprefix = 'simy';
day=1;
ep=2;
cellcountthresh = 5; % at least 5 cells fired as a candidate event
savedata = 1; % save data = 1; not save = 0
figopt = 1 % 1 = plot decoding result for each event; 0 = no figure generated
shuffleIterations = 0 % 1500 for standard decoding

figopt = 1 % 1 = plot decoding result for each event; 0 = no figure generated
shuffleIterations = 1500 % 1500 for standard decoding

preplay_decoding_CA1_singleday(animalprefix,day,ep,cellcountthresh, savedir, savedata, figopt, shuffleIterations)


pvals = replaytrajectory{day}{ep}.pvalue(:,1);
figure; histogram(pvals, 10)




%% Below, code that calls replay decoding for replayNet
clc
clear all
close all

%if ispc
%    userDir = 'C:/Users/Jordan/';
%elseif ismac
%    userDir = '/Users/jordanbreffle/';
%end
userDir = char(java.lang.System.getProperty('user.home'));

addpath(genpath([userDir, 'Box Sync/Code/Replay Project/Spiking network, original/sjlab_tools']))

simulationFileName = mfilename;
simulationFilePath = mfilename('fullpath');
startDateTime = clock; % Date and time the script was started
matlabVersion = version; % version of matlab used
systemVersion = ver; % Versions of matlab, toolboxes, Java, and OS
rngState = rng;       % the state of the rng 
tic
%%
%animalprefix_list = {'KL8'};          
animalprefix_list = {'2021_06_08_n1s1'};

%savedir = ('/Users/wenbotang/Bayesian_decoding/');% directory for saving results
savedir = [userDir, '/Box Sync/Data/Replay project/jb_simulation/v2_8_5/Bayesian Decode/'];% directory for saving results
runParamName = 'runParams_06_08_n1s1.mat'; %

eps = 2:2:2; 
day_list = 1;
cellcountthresh = 5; % at least 5 cells fired as a candidate event
savedata = 1; % save data = 1; not save = 0
figopt = 0 % 1 = plot decoding result for each event; 0 = no figure generated
shuffleIterations = 0 % 1500 for standard decoding

%savedata = 0; % save data = 1; not save = 0
%figopt = 2; % 1 = plot decoding result for each event; 0 = no figure generated, 2=make new fig for every event

%{
animalprefix_list = {'KL8'};
ep = 10; 
savedir0 = ('C:/Users/jtbre/Box Sync/Data/Replay project/js_SingleDayExpt/KL8ep10_Bayes_decode/');% directory for saving results
numReps = 50;
eps = repelem(ep, 1, numReps); 
epCounter = 1;
%}

%%
warning off
if ~exist(savedir)
    mkdir(savedir)
end
for theAnimal = 1:length(animalprefix_list)
    animalprefix = animalprefix_list{theAnimal};
    for day = day_list
        for ep = eps
            
            %savedir = sprintf('%sRep%s_', savedir0, num2str(epCounter)) % rm           
            
            disp(['Animal: ', animalprefix, ', Day: ', num2str(day), ', Episode: ', num2str(ep), ' Runtime: ', num2str(toc/60/60, 3), ' hours' ]);
            replay_decoding_CA1_singleday(animalprefix,day,ep,cellcountthresh,savedir,savedata, figopt, shuffleIterations)
            
            %epCounter = epCounter+1; % rm
            
        end
    end
end
warning on
runTime = toc
%%
save( strcat(savedir, runParamName), 'runTime', 'rngState', 'simulationFileName', ...
    'simulationFilePath', 'startDateTime', 'matlabVersion', 'systemVersion', ...
    'animalprefix_list', 'day_list', 'eps', 'cellcountthresh')

%rmpath(genpath('C:/Users/jtbre/Box Sync/Code/Replay Project/Spiking network, original/sjlab_tools'));
rmpath(genpath([userDir, 'Box Sync/Code/Replay Project/Spiking network, original/sjlab_tools']))

