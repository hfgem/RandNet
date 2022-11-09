%RandNet Project - random networks with global inhibition and their ability
%to produce sequences.
%
% randnet_PF is the same as randnet.m, but includes simulation of place
% fields for each network.
% 


%% Simulation options and save path

clear all

parameters.saveFlag = 0; % 1 to save simulation results
parameters.selectPath = 1; % 1 to select save destination, 0 to save in current dir
parameters.plotResults = 1; % 1 to plot basic simulation results

if parameters.saveFlag & parameters.selectPath
    save_path = uigetdir('/Users/hannahgermaine/Documents/PhD/','Select Save Folder'); %Have user input where they'd like the output stored
else
    save_path = [pwd, '/results/randnet_PF_Envs'];
end
addpath('functions')

[s,sim_git_hash_string] = system('git rev-parse HEAD');
parameters.sim_git_hash_string = sim_git_hash_string;

%% Initialize parameters

% Network structure parameters
parameters.n = 500; %number of neurons
parameters.clusters = 10; % Number of clusters in the network
parameters.mnc = 2; % mean number of clusters each neuron is a member of

% Time
parameters.t_max = 3; %maximum amount of time (s)
parameters.dt = 0.1*10^(-3); %timestep (s)

% Basic model parameters
% tau_E ~= 10 ms from direct data, DOI: 10.1126/science.aaf1836
parameters.tau_syn_E = 10*10^(-3); % Exc. synaptic decay time constant (s) PF19=50ms, HF18=10ms for figs 7-8 and longer for earlier figs
% tau_I ~= 1.2-8 ms from direct data, https://doi.org/10.1073/pnas.192233099
parameters.tau_syn_I = 3*10^(-3);  % Inh. synaptic decay time constant (s) PF19=5ms,  HF18=10ms for figs 7-8 and for earlier figs
parameters.tau_stdp = 5*10^(-3); %STDP time constant (s)                 
parameters.E_K = -80*10^(-3); %potassium reversal potential (V) %-75 or -80 mV
parameters.E_L = -70*10^(-3); %leak reversal potential (V) %-60 - -70 mV range
parameters.G_L = 25*10^(-9); %leak conductance (S) %10 - 30 nS range
parameters.C_m = 0.4*10^(-9); %total membrane capacitance (F) %Huge range from 0.1 - 100 pF
parameters.V_m_noise = 0.0*10^(-3); % 10^(-4); %magnitude of noise to use in membrane potential simulation (V)
parameters.V_th = -50*10^(-3); %threshold membrane potential (V)
parameters.V_reset = -70*10^(-3); %reset membrane potential (V)
parameters.V_syn_E = 0; %synaptic reversal potential (excitatory)
parameters.V_syn_I = -70*10^(-3); %synaptic reversal potential (inhibitory) %generally -70 pr -80 mV

% Recurrent connection strengths
parameters.del_G_syn_E_E = 750*10^(-12); %synaptic conductance step following spike (S)
parameters.del_G_syn_I_I = 0; %1.4*del_G_syn_E_E; %synaptic conductance step following spike (S)
parameters.del_G_syn_E_I = 500*10^(-12); %synaptic conductance step following spike (S)
parameters.del_G_syn_I_E = nan; %synaptic conductance step following spike (S)

% SRA parameters
parameters.del_G_sra = 330e-09; %spike rate adaptation conductance step following spike %ranges from 1-200 *10^(-9) (S)
parameters.tau_sra = 30*10^(-3); %spike rate adaptation time constant (s)

% STDP parameters
parameters.connectivity_gain = 0; %0.005; %amount to increase or decrease connectivity by with each spike (more at the range of 0.002-0.005)

% Input parameters:
% Poisson input
parameters.usePoisson = 1; % 1 to use poisson spike inputs, 0 for randn() input
parameters.rG = 1000; % input spiking rate, if using poisson inputs
parameters.W_gin = 750*10^-12; % increase in conductance, if using poisson inputs
parameters.cueSigma = 0; % temp value, to produce identical values
% Conductance input
parameters.G_std = -19*10^-9; % STD of the input conductance G_in, if using randn()
parameters.G_mean = 0* 10^-12; % mean of the input conductance G_in, if using randn()

% Network connection parameters
parameters.conn_prob = 0.08; %set a total desired connection probability
parameters.p_E = 0.75; %probability of an excitatory neuron
parameters.include_all = 2; % if a neuron is not in any cluster, take cluster membership from a highly connected neuron
parameters.global_inhib = 1; % if 1, I-cells are not clustered and have connection probability p_I
parameters.p_I = 0.5; % probability of an I cell connecting to any other cell

% Number of trials per net to run
parameters.nTrials = 1; % How many tests of different initializations to run
parameters.nNets = 1; % How many networks to run



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Temp, parameter tuning: %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%

PFsimFlag = 1;
PFscoreFlag = 0;
preplaySimFlag = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Combined Preplay and PFs
parameters.t_max = 12; %maximum amount of time (s)
parameters.G_L = 10*10^(-9); %leak conductance (S) %10 - 30 nS range
parameters.rG = 5000;
parameters.Win_var = (5e-12)^2;
parameters.tau_sra = 30*10^(-3); %spike rate adaptation time constant (s)
parameters.del_G_sra = 3.0e-012; %spike rate adaptation conductance step following spike %ranges from 1-200 *10^(-9) (S)
parameters.include_all = 3; % if a neuron is not in any cluster, take cluster membership from a highly connected neuron
parameters.clusters = 8; % Number of clusters in the network
parameters.mnc = 1.5; % mean number of clusters each neuron is a member of


parameters.Win_mean = 73 *10^-12;
parameters.IcueScale_PF = 1.4; % scales strength of I cell cue input, if ~=1 then Icells receive no spatial inputs
parameters.IcueScale = 1.01; % scales strength of I cell cue input, if ~=1 then Icells receive no spatial inputs
parameters.del_G_syn_E_E = 135*10^(-12); %synaptic conductance step following spike (S)
parameters.del_G_syn_E_I = 80*10^(-12); %synaptic conductance step following spike (S)
parameters.p_I = 0.25; % probability of an I cell connecting to any other cell
parameters.PFcontextFrac = 0.1; % scales E-cell's context cue input during PF trials, and 1-PFcontextFrac for spatial inputs


parameters.t_max = 60; %maximum amount of time (s)


parameters.del_G_syn_E_E = 150*10^(-12); %synaptic conductance step following spike (S)
parameters.del_G_syn_E_I = 90*10^(-12); %synaptic conductance step following spike (S)

% parameters.del_G_syn_E_E = 135*10^(-12); %synaptic conductance step following spike (S)
% parameters.del_G_syn_E_I = 75*10^(-12); %synaptic conductance step following spike (S)

parameters.clusterCorrs = 0; % old default
parameters.inhSRA = 1; % old default


parameters.clusterCorrs = 1
parameters.inputBiasSigma = 25

%{
parameters.Win_var = (8e-12)^2;
parameters.Win_mean = 73 *10^-12;
parameters.del_G_syn_E_E = 120*10^(-12); %synaptic conductance step following spike (S)
parameters.del_G_syn_E_I = 150*10^(-12); %synaptic conductance step following spike (S)
%}


%{
% Uncomment for example network visualization
% Then run randnet_plotting.m network section
PFsimFlag = 0;
PFscoreFlag = 0;
preplaySimFlag = 0;
parameters.clusters = 8
parameters.mnc = 1.5
%}


%% Attempt at minimal changes for IS regime
%{
parameters.inhSRA = 0; % no SRA for I cells

parameters.del_G_syn_E_E = 0*180*10^(-12); %synaptic conductance step following spike (S)
parameters.del_G_syn_E_I = 0*350*10^(-12); %synaptic conductance step following spike (S)

parameters.Win_mean = 71 *10^-12;
parameters.IcueScale_PF = 0.75;
parameters.IcueScale = 0.75; % scales strength of I cell cue input, if ~=1 then Icells receive no spatial inputs

PFsimFlag = 0;
PFscoreFlag = 0;
preplaySimFlag = 1 ;
parameters.t_max = 6; %maximum amount of time (s)
%}

%% Add in cluster correlations to location input

%{
parameters.clusterCorrs = 1;
parameters.inputBiasSigma = 25; 
parameters.inhSRA = 0; 

parameters.del_G_syn_E_E = 135*10^(-12); %synaptic conductance step following spike (S)
parameters.del_G_syn_E_I = 90*10^(-12); %synaptic conductance step following spike (S)

parameters.clusters = 10; % Number of clusters in the network
parameters.mnc = 1.5; % mean number of clusters each neuron is a member of


parameters.Win_mean = 72 *10^-12;
parameters.Win_var = (5e-12)^2;

PFsimFlag = 1;
PFscoreFlag = 0;
preplaySimFlag = 1;


parameters.mnc = 1.5; % mean number of clusters each neuron is a member of
parameters.del_G_syn_E_E = 220*10^(-12); %synaptic conductance step following spike (S)
parameters.del_G_syn_E_I = 400*10^(-12); %synaptic conductance step following spike (S)
parameters.t_max = 6;


parameters.IcueScale_PF = 1.0; % scales strength of I cell cue input, if ~=1 then Icells receive no spatial inputs
parameters.IcueScale =  0.75; % scales strength of I cell cue input, if ~=1 then Icells receive no spatial inputs
%}

%% Attempt at IS regime, rate of 500, no SRA for Inh cells

% parameters.clusters = 10; % Number of clusters in the network
% parameters.mnc = 1.5; % mean number of clusters each neuron is a member of
% parameters.del_G_sra = 3.0e-012; %spike rate adaptation conductance step following spike %ranges from 1-200 *10^(-9) (S)
% parameters.p_I = 0.5; % probability of an I cell connecting to any other cell

%{
parameters.inhSRA = 0; 
%} 
%{
parameters.del_G_syn_E_E = 400*10^(-12); %synaptic conductance step following spike (S)
parameters.del_G_syn_E_I = 500*10^(-12); %synaptic conductance step following spike (S)
parameters.Win_mean = 600 *10^-12;
parameters.rG = 500;
parameters.IcueScale = 1.01;


parameters.Win_mean = 590 *10^-12;
parameters.del_G_syn_E_E = 470 *10^(-12); %synaptic conductance step following spike (S)
parameters.del_G_syn_E_I = 525 *10^(-12); %synaptic conductance step following spike (S)
parameters.clusters = 6; % Number of clusters in the network
parameters.mnc = 1.25; % mean number of clusters each neuron is a member of


parameters.Win_mean = 585 *10^-12;
parameters.del_G_syn_E_E = 600 *10^(-12); %synaptic conductance step following spike (S)
parameters.del_G_syn_E_I = parameters.del_G_syn_E_E; %synaptic conductance step following spike (S) E to I
parameters.del_G_syn_I_E = 2.25 *parameters.del_G_syn_E_I; %synaptic conductance step following spike (S) I to E

parameters.clusters = 8; % Number of clusters in the network
parameters.mnc = 1.5; % mean number of clusters each neuron is a member of
%}

%{
parameters.Win_mean = 585 *10^-12;
parameters.del_G_syn_E_E = 620 *10^(-12); %synaptic conductance step following spike (S)
parameters.del_G_syn_E_I = 630 *10^(-12); %synaptic conductance step following spike (S) E to I
parameters.del_G_syn_I_E = 2.4* 600 *10^(-12); %synaptic conductance step following spike (S) I to E
parameters.del_G_sra = 100.0e-012;
%}


%% Attempt at IS regime, rate of 1000, no SRA for Inh cells

%{
parameters.inhSRA = 0; 
parameters.t_max = 6;

parameters.rG = 1000;
parameters.Win_mean = 315 *10^-12;

parameters.del_G_syn_E_E = 350*10^(-12); %synaptic conductance step following spike (S)
parameters.del_G_syn_E_I = 350*10^(-12); %synaptic conductance step following spike (S)
parameters.del_G_syn_I_E = 3.0* parameters.del_G_syn_E_I; %synaptic conductance step following spike (S) I to E
parameters.IcueScale = 0.75;

parameters.Win_mean = 325 *10^-12;
parameters.IcueScale_PF = 0.90; % scales strength of I cell cue input, if ~=1 then Icells receive no spatial inputs
%}
%{
parameters.del_G_syn_E_E = 350*10^(-12); %synaptic conductance step following spike (S)
parameters.del_G_syn_E_I = 2.5 *350*10^(-12); %synaptic conductance step following spike (S)
parameters.del_G_syn_I_E = 1.0 *350*10^(-12); %synaptic conductance step following spike (S) I to E
parameters.IcueScale = 0.5;
%}
%{
parameters.del_G_syn_E_E = 400*10^(-12); %synaptic conductance step following spike (S)
parameters.del_G_syn_E_I = 2.75 *350*10^(-12); %synaptic conductance step following spike (S)
parameters.del_G_syn_I_E = 0.75 *350*10^(-12); %synaptic conductance step following spike (S) I to E
parameters.IcueScale = 0.75;
%}


%% Lower GL, rG=5000, no SRA for Inh cells

%{
parameters.inhSRA = 0; 

parameters.G_L = 5*10^(-9); %leak conductance (S) %10 - 30 nS range

parameters.Win_mean = 37 *10^-12;

parameters.del_G_syn_E_E = 150*10^(-12); %synaptic conductance step following spike (S)
parameters.del_G_syn_E_I = 150*10^(-12); %synaptic conductance step following spike (S)
parameters.del_G_syn_I_E = 1.5* parameters.del_G_syn_E_I; %synaptic conductance step following spike (S) I to E
parameters.IcueScale = 0.95;

parameters.Win_mean = 37 *10^-12;
parameters.del_G_syn_E_E = 0.5*parameters.del_G_syn_E_E; %synaptic conductance step following spike (S)
parameters.IcueScale_PF = 1.4;


PFsimFlag = 1;
PFscoreFlag = 0;
preplaySimFlag = 0;
parameters.t_max = 1.5; %maximum amount of time (s)
%}

%% Parameters for sequence analysis

% Analysis parameters for PBE detection
parameters.PBE_min_Hz = 0.5; % minimum population mean rate during PBE
parameters.PBE_zscore = 1.0; % minimum stds above mean rate to detect PBE
parameters.PBE_min_dur = 30 * (1/1000); % minimum duration of a PBE
parameters.PBE_window =  15 * (1/1000) *(1/parameters.dt); % width of gaussian kernel used to calculate mean pop activity rate
parameters.PBE_max_combine = 10 * (1/1000); % Combine adjacent PBEs separaeted by less than this duration

% parameters.PBE_zscore = 1.0; % minimum stds above mean rate to detect PBE
% parameters.returnToMean = 1; % if 1, terminate events at leave/return to mean; else, leave/return to PBE_zscore

%% __set/update Dependent Parameters__ %%

parameters = set_depedent_parameters(parameters);

%Save to computer
if parameters.saveFlag == 1
    save(strcat(save_path,'/parameters.mat'),'parameters')
end


%% Place field simulation parameters

% Copy preplay parameters as defaults
pfsim = parameters; 
pfsim.PFscoreFlag = PFscoreFlag;

% Add PF specific parameters and PF specific values
pfsim.envIDs = [1, 2, 3, 4] ; % odd IDs are rightward traversals, even IDs are rightward traversals of the ID-1 track
pfsim.nEnvironments = numel(pfsim.envIDs);
pfsim.nTrials = 5;
pfsim.t_max = 2;
pfsim.t = 0:parameters.dt:pfsim.t_max;

pfsim.xPos = 0:parameters.dt/pfsim.t_max:1;
pfsim.trackWidth = 1; % maximum value
pfsim.spatialBin = 2/100; % 2 cm bins for localizing position
pfsim.linFieldGaussSD = 0.04;% Standard deviation of PF gaussian kernel
pfsim.winScale = 5; % window is 5x the STDev
pfsim.gridxvals = pfsim.spatialBin:pfsim.spatialBin:pfsim.trackWidth;          % Grid-points in the x-direction

% PF Gaussian fit parameters
pfsim.gaussFOLower = [10, 0, 0]; % [peak amplitude, position of peak on track, standard deviation of peak]
pfsim.gaussFOUpper = [30, max(pfsim.gridxvals)*100, sqrt(max(pfsim.gridxvals)*100)]; 
pfsim.peakTarget = 15; % Hz, Target peak rate for linfieldsScore

pfsim.minPeakRate = 3; % minimum PF peak rate to consider a cell a place cell

% set depedendent parameters for PF sim
pfsim = set_depedent_parameters(pfsim);

%Save to computer
if parameters.saveFlag == 1
    save(strcat(save_path,'/pfsim.mat'),'pfsim')
end

%% Create Networks and Check Spike Progression
%Runs through a series of different random number generator seeds to change
%the network connectivity and setup, and then automatically outputs
%sequence data to a folder in save_path

%If uploading a parameter file, uncomment the next line
% load(strcat(save_path,'/parameters.mat'))

if parameters.saveFlag & ~isfolder(save_path)
    mkdir(save_path);
end
    
if parameters.saveFlag
    save(strcat(save_path,'/parameters.mat'),'parameters'); 
end

%____________________________________________
%___Run Simulations for Different Networks___
%____________________________________________

for ithNet = 1:parameters.nNets
    
    %CREATE NETWORK SAVE PATH
    net_save_path = strcat(save_path,'/network_',string(ithNet));
    if parameters.saveFlag & ~isfolder(net_save_path)
        mkdir(net_save_path);
    end
    
    % Generete ith network structure
    network = create_clusters(pfsim, 'seed', ithNet, 'include_all', parameters.include_all, 'global_inhib', parameters.global_inhib);
    if parameters.saveFlag
        save(strcat(net_save_path,'/network.mat'),'network');
    end
    
    
    %% Place field simulation:
    startTime = now;
    tic;
    if PFsimFlag
        opV = zeros(parameters.n, numel(pfsim.t), pfsim.nEnvironments, pfsim.nTrials); % Voltage from all sims
        opS = zeros(parameters.n, numel(pfsim.t), pfsim.nEnvironments, pfsim.nTrials); % Spikes from all sims
        
        for ithEnv = 1:pfsim.nEnvironments
            
            for ithTrial = 1:pfsim.nTrials
                
                G_in_PFs = zeros(parameters.n, numel(pfsim.t));
                %rng(ithTrial)
                % G_in_PFs(:,1,ithEnv,ithTrial) = 1/10* dI(:,ithEnv) * 2*parameters.rGmax * parameters.tau_syn_E + sqrt(1/2*parameters.tau_syn_E*dI(:,ithEnv).^2*2*parameters.rGmax).*randn(parameters.n, 1) ; 
                G_in_PFs(:,1) = zeros(parameters.n, 1); 
                for i = 2:numel(pfsim.t)
                        % Exponential decay from last time step
                        G_in_PFs(:,i) = G_in_PFs(:,i-1)*exp(-parameters.dt/parameters.tau_syn_E);
                        % Add spikes from each input source
                        G_in_PFs(:,i) = G_in_PFs(:,i) + ...
                                network.spatialInput{1}(:,ithEnv) .* (1-parameters.PFcontextFrac).*[rand(parameters.n, 1) < (parameters.dt* (parameters.rG * (i/numel(pfsim.t)) ))] + ...
                                network.spatialInput{2}(:,ithEnv) .* (1-parameters.PFcontextFrac).* [rand(parameters.n, 1) < (parameters.dt* ( parameters.rG * ((numel(pfsim.t)-i)/numel(pfsim.t)) ) )] + ...
                                network.contextInput(:,ithEnv) .* [parameters.PFcontextFrac.*ismember(network.all_indices, network.E_indices)]' .* [rand(parameters.n, 1) < (parameters.dt*parameters.rG)] + ...
                                network.contextInput(:,ithEnv) .* [parameters.IcueScale_PF.*ismember(network.all_indices, network.I_indices)]' .* [rand(parameters.n, 1) < (parameters.dt*parameters.rG)]   ;
                end
                
                % Set up for simulation
                V_m = zeros(parameters.n,numel(pfsim.t)); %membrane potential for each neuron at each timestep
                V_m(:,1) = -60e-3 + 5e-3*randn([parameters.n,1]); %set all neurons to baseline reset membrane potential with added noise
                pfsim.G_in = G_in_PFs; 
                trialSeed = randi(10^6);

                % PF Simulation
                [V_m, G_sra, G_syn_E_E, G_syn_I_E, G_syn_E_I, G_syn_I_I, conns] = randnet_calculator(pfsim, trialSeed, network, V_m);
                opV(:,:,ithEnv,ithTrial) = V_m;
                opS(:,:,ithEnv,ithTrial) = logical( V_m>parameters.V_th);
                
            end
            

        end
        % keyboard

        % Calculate Place fields and PF-objective score and plot PF sim figs
        allScores = zeros(1, pfsim.nEnvironments);
        % allPFpeaksSeq = [];
        % allLinfields = [];
        
        [linfields, PFpeaksSequence] = calculate_linfields(opS, pfsim, pfsim, network, true);
        
        for ithEnv = 1:pfsim.nEnvironments
            %{
            trialSpikes = opS(:,:,ithEnv,:);
            [linfields, PFpeaksSequence] = calculate_linfields(trialSpikes, pfsim, pfsim, network, true);
            allPFpeaksSeq(:,ithEnv) = PFpeaksSequence;
            allLinfields{ithEnv} = linfields;
            %}
            
            %{
            figure; plot(pfsim.t, V_m(network.E_indices(1:3),:)); ylabel('Vm (V)'); xlabel('Time (s)'); 
            figure; plot(pfsim.t, V_m(network.I_indices(1:3),:)); ylabel('Vm (V)'); xlabel('Time (s)'); 
            figure; plot(pfsim.t, G_in_PFs(1:4,:,1,1)); ylabel('G in (S)'); xlabel('Time (s)');      
            figure; plotSpikeRaster( logical( [ opS(network.E_indices,:,1,1), opS(network.I_indices,:,1,1) ]) ...
            'TimePerBin', parameters.dt, 'PlotType', 'scatter'); % 
            ylabel('Cell');
            figure; plotSpikeRaster( logical(opS(network.I_indices,:,1,1)), ...
            'TimePerBin', parameters.dt, 'PlotType', 'scatter'); % 
            ylabel('Cell');
            %}
            if pfsim.PFscoreFlag % TODO: adapt to {ithEnv} in {traj} dimension
                allScores(ithEnv) = calculate_linfieldsScore(linfields, pfsim, pfsim, network)
                disp('Has calculate_linfieldsScore been fixed for multiple environments?')
            end
            disp(['Env: ', num2str(ithEnv), ', Score: ', num2str(allScores(ithEnv))])
            
            %figure; plotSpikeRaster( logical( [ opS(network.E_indices,:,ithEnv,1); opS(network.I_indices,:,1,1) ]), 'TimePerBin', parameters.dt, 'PlotType', 'scatter');
            rpermIcells = randperm(numel(network.I_indices));
            figure; plotSpikeRaster( logical( [ opS(network.E_indices(PFpeaksSequence{ithEnv}),:,ithEnv,1); opS(network.I_indices(rpermIcells),:,ithEnv,1) ]), 'TimePerBin', parameters.dt, 'PlotType', 'scatter');

        end
        score = mean(allScores);
        disp(['Mean score: ', num2str(score)])

        figure; histogram( sum(opS, [2:4])./parameters.t_max./parameters.nTrials./pfsim.nEnvironments, 50 )
        xlabel('Mean rate (Hz, all PF trials)'); ylabel('All cells')
    
        plotPFs(linfields, PFpeaksSequence, network, pfsim)
        
    end
    PFruntime = toc
    
        
    %% Preplay simulation:
    tic
    if preplaySimFlag
        V_m_var = struct;
        G_var = struct;
        I_var = struct;
        network_spike_sequences = struct;
        network_cluster_sequences = struct; 

        for ithTest = 1:parameters.nTrials       

            %Create input conductance variable
            rng(ithTest)
            preplayContext = 1; % Which environments context cue to use for preplay sim
            if parameters.usePoisson
                G_in = zeros(parameters.n, parameters.t_steps+1);
                G_in(:,1) = 1/2* parameters.Win_mean * 2*parameters.rG * parameters.tau_syn_E + sqrt(1/2*parameters.tau_syn_E*parameters.Win_mean.^2*2*parameters.rG).*randn(parameters.n, 1) ; 
                for k = 2:(parameters.t_steps+1)
                    G_in(:,k) = G_in(:,k-1)*exp(-parameters.dt/parameters.tau_syn_E);

                    % G_in(:,k) = G_in(:,k) + network.contextInput .* [rand(parameters.n, 1) < (parameters.dt*parameters.rG)];
                    G_in(:,k) = G_in(:,k) + [network.contextInput(:,preplayContext) .* [1     .*              ismember(network.all_indices, network.E_indices)]' .* [rand(parameters.n, 1) < (parameters.dt*parameters.rG)] + ...
                                             network.contextInput(:,preplayContext) .* [parameters.IcueScale.*ismember(network.all_indices, network.I_indices)]'  .* [rand(parameters.n, 1) < (parameters.dt*parameters.rG)]] ;
                end
            else
                G_in = (parameters.G_std*randn(parameters.n,parameters.t_steps+1))+parameters.G_mean;
                G_in(G_in<0) = 0;
            end
            parameters.('G_in') = G_in;

            %Create Storage Variables
            V_m = zeros(parameters.n,parameters.t_steps+1); %membrane potential for each neuron at each timestep
            V_m(:,1) = -60e-3 + 5e-3*randn([parameters.n,1]); %set all neurons to baseline reset membrane potential with added noise
            %indSpike = randperm(numel(network.E_indices)); indSpike = indSpike(1:20);
            %V_m(network.E_indices(indSpike),1) = parameters.V_th-0.01e-3;
            %V_m(:,1) = parameters.V_reset + randn([parameters.n,1])*(10^(-3))*sqrt(parameters.dt); %set all neurons to baseline reset membrane potential with added noise

            seed = ithTest;

            %Run model
            [V_m, G_sra, G_syn_E_E, G_syn_I_E, G_syn_E_I, G_syn_I_I, conns] = ...
                    randnet_calculator(parameters, seed, network, V_m);
            V_m_var(ithTest).V_m = V_m;
            G_var(ithTest).G_in = G_in;
            G_var(ithTest).G_sra = G_sra;
            G_var(ithTest).G_syn_I_E = G_syn_I_E;
            G_var(ithTest).G_syn_E_E = G_syn_E_E;
            G_var(ithTest).G_syn_I_I = G_syn_I_I;
            G_var(ithTest).G_syn_E_I = G_syn_E_I;

            %{
            [V_m, conns] = ...
                    randnet_calculator_memOpt(parameters, seed, network, V_m);
            V_m_var(ithTest).V_m = V_m;
            %}

            [trialResults] = detect_PBE( V_m(network.E_indices,:)>= parameters.V_th, parameters);
            if ithTest == 1 % append trialResults struct to network results struct
                network_spike_sequences = trialResults;
            else
                network_spike_sequences = [network_spike_sequences, trialResults]; 
            end

            if parameters.plotResults
                plot_randnet_results(parameters, network, V_m, G_in, network_spike_sequences, ithTest, net_save_path)
            end


        end % trial loop
    end
    Preplayruntime = toc

    
    %SAVE NETWORK DATA
    if parameters.saveFlag
        save(strcat(net_save_path,'/V_m_var.mat'),'V_m_var','-v7.3')
        save(strcat(net_save_path,'/G_var.mat'),'G_var','-v7.3')
        save(strcat(net_save_path,'/I_var.mat'),'I_var','-v7.3')
        save(strcat(net_save_path,'/network_spike_sequences.mat'),'network_spike_sequences','-v7.3')
        save(strcat(net_save_path,'/network_cluster_sequences.mat'),'network_cluster_sequences','-v7.3')
    end
    
end % Network loop

