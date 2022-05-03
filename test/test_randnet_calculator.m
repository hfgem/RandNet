% Tests that randnet_calculator and randnet_calculator_memOpt produce 
% identical results


%% Initialize parameters

% Network structure parameters
parameters.n = 500; %number of neurons
parameters.clusters = 10; % Number of clusters in the network
parameters.mnc = 2; % mean number of clusters each neuron is a member of

% Time
parameters.t_max = 1; %maximum amount of time (s)
parameters.dt = 0.1*10^(-3); %timestep (s)

% Basic model parameters
parameters.tau_syn_E = 10*10^(-3); % Exc. synaptic decay time constant (s) PF19=50ms, HF18=10ms for figs 7-8 and longer for earlier figs
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

% Network connection parameters
parameters.conn_prob = 0.08; %set a total desired connection probability
parameters.p_E = 0.75; %probability of an excitatory neuron
parameters.include_all = 1; % if a neuron is not in any cluster, take cluster membership from a highly connected neuron
parameters.global_inhib = 1; % if 1, I-cells are not clustered and have connection probability p_I
parameters.p_I = 0.5; % probability of an I cell connecting to any other cell

% Number of trials per net to run
parameters.nTrials = 1; % How many tests of different initializations to run
parameters.nNets = 1; % How many networks to run


%% __set/update Dependent Parameters__ %%
parameters = set_depedent_parameters(parameters);


%% Generete ith network structure
ithNet = 1;
network = create_clusters(parameters, 'seed', ithNet, 'include_all', parameters.include_all, 'global_inhib', parameters.global_inhib);


%% Create input conductance variable
if parameters.usePoisson
    G_in = single(zeros(parameters.n, parameters.t_steps+1));
    %G_in = zeros(parameters.n, parameters.t_steps+1);
    for k = 2:(parameters.t_steps+1)
        G_in(:,k) = G_in(:,k-1)*exp(-parameters.dt/parameters.tau_syn_E);
        G_in(:,k) = G_in(:,k) + parameters.W_gin * [rand(parameters.n, 1) < (parameters.dt*parameters.rG)];
    end
else
    G_in = (parameters.G_std*randn(parameters.n,parameters.t_steps+1))+parameters.G_mean;
    G_in(G_in<0) = 0;
end
parameters.('G_in') = G_in;


%% Create Storage Variables
V_m = zeros(parameters.n,parameters.t_steps+1); %membrane potential for each neuron at each timestep
V_m(:,1) = parameters.V_reset + randn([parameters.n,1])*(10^(-3))*sqrt(parameters.dt); %set all neurons to baseline reset membrane potential with added noise
seed = 1;


%% Run model

[V_m1] = randnet_calculator(parameters, seed, network, V_m);
spikes1 = logical(V_m1>=parameters.V_th);

[spikes2] = randnet_calculator_memOpt(parameters, seed, network, V_m(:,1) );
        

%% Test

assert( isequal(spikes1, spikes2), 'Not Equal' )
isequal(spikes1, spikes2)















