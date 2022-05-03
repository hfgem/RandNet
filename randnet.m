%RandNet Project - random networks with global inhibition and their ability
%to produce sequences

%% Simulation options and save path

clear all

parameters.saveFlag = 0; % 1 to save simulation results
parameters.selectPath = 0; % 1 to select save destination, 0 to save in current dir
parameters.plotResults = 1; % 1 to plot basic simulation results

if parameters.saveFlag & parameters.selectPath
    save_path = uigetdir('/Users/hannahgermaine/Documents/PhD/','Select Save Folder'); %Have user input where they'd like the output stored
else
    save_path = [pwd, '/results/randnet'];
end
addpath('functions')

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
parameters.tau_syn_I = 2*10^(-3);  % Inh. synaptic decay time constant (s) PF19=5ms,  HF18=10ms for figs 7-8 and for earlier figs
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
parameters.del_G_syn_I_E = 500*10^(-12); %synaptic conductance step following spike (S)

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
parameters.include_all = 1; % if a neuron is not in any cluster, take cluster membership from a highly connected neuron
parameters.global_inhib = 1; % if 1, I-cells are not clustered and have connection probability p_I
parameters.p_I = 0.5; % probability of an I cell connecting to any other cell

% Number of trials per net to run
parameters.nTrials = 1; % How many tests of different initializations to run
parameters.nNets = 1; % How many networks to run


%{
% Alternative working parameter set
parameters.del_G_syn_E_E = 675*10^(-12); %synaptic conductance step following spike (S)
parameters.del_G_syn_I_I = 0; %1.4*del_G_syn_E_E; %synaptic conductance step following spike (S)
parameters.del_G_syn_E_I = 450*10^(-12); %synaptic conductance step following spike (S)
parameters.del_G_syn_I_E = 450*10^(-12); %synaptic conductance step following spike (S)
parameters.del_G_sra = 100e-09; % or 50nA spike rate adaptation conductance step following spike %ranges from 1-200 *10^(-9) (S)
%}

%{
% Attempt with lower dSRA
parameters.del_G_syn_E_E = 575*10^(-12); %synaptic conductance step following spike (S)
parameters.del_G_syn_I_I = 0; %1.4*del_G_syn_E_E; %synaptic conductance step following spike (S)
parameters.del_G_syn_E_I = 325*10^(-12); %synaptic conductance step following spike (S)
parameters.del_G_syn_I_E = 325*10^(-12); %synaptic conductance step following spike (S)
parameters.del_G_sra = 10e-09; %spike rate adaptation conductance step following spike %ranges from 1-200 *10^(-9) (S)
%}

%{
% Attempt with lower rG
parameters.rG = 250; % input spiking rate, if using poisson inputs
parameters.W_gin = 2300*10^-12; % increase in conductance, if using poisson inputs
parameters.del_G_syn_E_E = 1200*10^(-12); %synaptic conductance step following spike (S)
parameters.del_G_syn_I_I = 0; %1.4*del_G_syn_E_E; %synaptic conductance step following spike (S)
parameters.del_G_syn_E_I = 550*10^(-12); %synaptic conductance step following spike (S)
parameters.del_G_syn_I_E = 550*10^(-12); %synaptic conductance step following spike (S)
%}


%% Parameters for sequence analysis

%{
parameters.E_events_only = 1; % if 1, only consider E-cells for detect_events
parameters.IEI = 0.02; %inter-event-interval (s) the elapsed time between spikes to count separate events
parameters.bin_width = 5*10^(-3); %5 ms bin

%TEST 1: The number of neurons participating in a sequence must pass a threshold:
parameters.event_cutoff = 0.10; %0.25; %fraction of neurons that have to be involved to constitute a successful event

%TEST 2: The firing rate must fall within a realistic range
parameters.min_avg_fr = 0.01;
parameters.max_avg_fr = 3.0;

% TEST 3: The sequence(s) of firing is(are) within reasonable lengths
parameters.min_avg_length = 0.01;
parameters.max_avg_length = 0.5;
%}

% Analysis parameters for PBE detection
parameters.PBE_min_Hz = 0.5; % minimum population mean rate during PBE
parameters.PBE_zscore = 1.0; % minimum stds above mean rate to detect PBE
parameters.PBE_min_dur = 30 * (1/1000); % minimum duration of a PBE
parameters.PBE_window =  10 * (1/1000) *(1/parameters.dt); % width of gaussian kernel used to calculate mean pop activity rate
parameters.PBE_max_combine = 10 * (1/1000); % Combine adjacent PBEs separaeted by less than this duration


%% __set/update Dependent Parameters__ %%
parameters = set_depedent_parameters(parameters);


%Save to computer
if parameters.saveFlag == 1
    save(strcat(save_path,'/parameters.mat'),'parameters')
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
    
    % Create ith networks savepath
    net_save_path = strcat(save_path,'/network_',string(ithNet));
    if parameters.saveFlag & ~isfolder(net_save_path)
        mkdir(net_save_path);
    end
    
    % Generete ith network structure
    network = create_clusters(parameters, 'seed', ithNet, 'include_all', parameters.include_all, 'global_inhib', parameters.global_inhib);
    if parameters.saveFlag
        save(strcat(net_save_path,'/network.mat'),'network');
    end
    
    
    %% Place field simulation:
    
    G_in_PFs = zeros(parameters.n, numel(sim.t), sim.nEnvironments_PFs, sim.nTrials_PFs);
    for ithTrial = 1:sim.nTrials
        G_in_PFs(:,1,ithEnv,ithTrial) = 1/10* dI(:,ithEnv) * 2*parameters.rGmax * parameters.tau_syn_E + sqrt(1/2*parameters.tau_syn_E*dI(:,ithEnv).^2*2*parameters.rGmax).*randn(parameters.n, 1) ; % ... and match variance

        for i = 2:numel(sim.t)
                % Exponential decay from last time step
                G_in_PFs(:,i,ithEnv,ithTrial) = G_in_PFs(:,i-1,ithEnv,ithTrial)*exp(-parameters.dt/parameters.tau_syn_E);

                G_in_PFs(:,i,ithEnv,ithTrial) = G_in_PFs(:,i,ithEnv,ithTrial) + ...
                        Wsig1 .* [rand(parameters.n, 1) < (parameters.dt* (parameters.rG * (i/numel(sim.t)) ))] + ...
                        Wsig2 .* [rand(parameters.n, 1) < (parameters.dt* ( parameters.rG * ((numel(sim.t)-i)/numel(sim.t)) ) )] + ...
                        WsigContext(:,ithEnv) .* [rand(parameters.n, 1) < (parameters.dt*parameters.rG)]  ;
            end
    end
    opV = zeros(parameters.n, numel(sim.t), sim.nEnvironments, sim.nTrials); % Voltage from all sims
    opS = zeros(parameters.n, numel(sim.t), sim.nEnvironments, sim.nTrials); % Spikes from all sims
    for ithEnv = 1:sim.nEnvironments
        for i = 1:sim.nTrials_PFs

            %Create Storage Variables
            I_syn = zeros(parameters.n,numel(sim.t)); %synaptic current emitted by each neuron at each timestep (A)
            %synaptic conductances for each neuron at each timestep
            V_m = zeros(parameters.n,numel(sim.t)); %membrane potential for each neuron at each timestep
            V_m(:,1) = -60e-3 + 1e-3*randn([parameters.n,1]); %set all neurons to baseline reset membrane potential with added noise
            parameters.G_in = G_in_PFs(:,:,ithEnv,ithTrial); 
            trialSeed = randi(10^6);
            
            [V_m, G_sra, G_syn_E_E, G_syn_I_E, G_syn_E_I, G_syn_I_I, conns] = randnet_calculator(parameters, trialSeed, network, V_m);
            opV(:,:,ithEnv,i) = V_m;
            opS(:,:,ithEnv,i) = V_m>parameters.V_th;
        end
    end
    
    % Calculate PF objective score
    for i = 1:sim.nEnvironments
        [allScores(i), linfields] = PFopt_ObjFun(opS(:,:,i,:), parameters, sim, track, network);
        disp(['Env: ', num2str(i), ', Score: ', num2str(allScores(i))])
    end
    score = mean(allScores);
    disp(['Mean score: ', num2str(score)])
    PF_sameNet_plottingv1_0

    % Extract place field order from linfields struct
    op = [];
    for epoch = 1:sim.nEnvironments
        opTemp = [];
        for ithCell = 1:parameters.n
            PF = linfields{day}{epoch}{tetrode}{ithCell}{tr}(:,5);
            %if sum(PF>0)
                opTemp= [opTemp;PF'];
            %end
        end
        op(epoch,:,:) = opTemp(network.E_indices,:);
    end
    row_all_zeros1 = find(all( squeeze(op(1,:,:))==0, 2)) ;
    row_n_all_zeros1 = find(~all( squeeze(op(1,:,:))==0, 2)) ;

    [peakRate,peakRateLocation] = max(squeeze(op(1,row_n_all_zeros1,:)), [], 2);
    [B,sortedCellIndsbyPeakRateLocation] = sort(peakRateLocation, 'descend');
    PFpeaksSequence = [row_n_all_zeros1(sortedCellIndsbyPeakRateLocation); row_all_zeros1];
    rateDenom1 = squeeze(max(op(1,PFpeaksSequence,:), [], 3))';
    minPeakRate = 3;
    PFpeaksSequence(rateDenom1<minPeakRate) = nan;

        
    %% Preplay simulation
    %Run through every cluster initialization and store relevant data and
    %calculations
    V_m_var = struct;
    G_var = struct;
    I_var = struct;
    network_spike_sequences = struct;
    
    for ithTest = 1:parameters.nTrials       
               
        % Create input conductance variable
        if parameters.usePoisson
            G_in = zeros(parameters.n, parameters.t_steps+1);
            for k = 2:(parameters.t_steps+1)
                G_in(:,k) = G_in(:,k-1)*exp(-parameters.dt/parameters.tau_syn_E);
                G_in(:,k) = G_in(:,k) + network.contextInput .* [rand(parameters.n, 1) < (parameters.dt*parameters.rG)];
            end
        else
            G_in = (parameters.G_std*randn(parameters.n,parameters.t_steps+1))+parameters.G_mean;
            G_in(G_in<0) = 0;
        end
        parameters.('G_in') = G_in;

        %Create Storage Variables
        V_m = zeros(parameters.n,parameters.t_steps+1); %membrane potential for each neuron at each timestep
        V_m(:,1) = parameters.V_reset + randn([parameters.n,1])*(10^(-3))*sqrt(parameters.dt); %set all neurons to baseline reset membrane potential with added noise
        
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

        % Detect PBEs from E-cell spikes
        [trialResults] = detect_PBE( V_m(network.E_indices,:)>= parameters.V_th, parameters);
        
        % Append trialResults struct to network results struct
        if ithTest == 1 
            network_spike_sequences = trialResults;
        else
            network_spike_sequences = [network_spike_sequences, trialResults]; 
        end
        
        if parameters.plotResults
            plot_randnet_results(parameters, network, V_m, G_in, network_spike_sequences, ithTest, net_save_path)
        end
        
        
    end % trial loop
    
    
    %SAVE NETWORK DATA
    if parameters.saveFlag
        save(strcat(net_save_path,'/V_m_var.mat'),'V_m_var','-v7.3')
        save(strcat(net_save_path,'/G_var.mat'),'G_var','-v7.3')
        save(strcat(net_save_path,'/I_var.mat'),'I_var','-v7.3')
        save(strcat(net_save_path,'/network_spike_sequences.mat'),'network_spike_sequences','-v7.3')
    end
    
end % Network loop

