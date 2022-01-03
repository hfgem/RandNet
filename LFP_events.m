%LFP Test - using the LFP to indicate a "SWR event" and look at the spikes
%there.
% Use a long simulation and search across the interval for areas of high
% fluctuation. Calculate LFP based on Mazzoni et al. 2015.

%% Save Path

clear all
save_path = uigetdir('/Users/hannahgermaine/Documents/PhD/','Select Save Folder'); %Have user input where they'd like the output stored

%% Initialization

%____________________________
%___Independent Parameters___
%____________________________

%Neuron and cluster counts
n = 200; %number of neurons
clusters = round(n/20); %number of clusters of neurons (for small n round(n/5), for large n round(n/20)) 

%Interaction constants
t_max = 10; %maximum amount of time (s)
dt = 0.1*10^(-3); %timestep (s)
tau_syn_E = 10*10^(-3); %AMPA/NMDA synaptic decay time constant (s)
tau_syn_I = 5*10^(-3); %GABA synaptic decay time constant (s)
tau_stdp = 5*10^(-3); %STDP time constant (s)
E_K = -80*10^(-3); %potassium reversal potential (V) %-75 or -80 mV
E_L = -70*10^(-3); %leak reversal potential (V) %-60 - -70 mV range
G_L = 25*10^(-9); %leak conductance (S) %10 - 30 nS range
C_m = 0.4*10^(-9); %total membrane capacitance (F) %Huge range from 0.1 - 100 pF
V_m_noise = 10^(-4); %magnitude of noise to use in membrane potential simulation (V)
V_th = -50*10^(-3); %threshold membrane potential (V)
V_reset = -70*10^(-3); %reset membrane potential (V)
V_syn_E = 0; %synaptic reversal potential (excitatory)
V_syn_I = -70*10^(-3); %synaptic reversal potential (inhibitory) %generally -70 pr -80 mV
%______Split del_G_syn______
% del_G_syn_E = 2*10^(-9); %synaptic conductance step following spike (S)
% del_G_syn_I = 8*10^(-9); %synaptic conductance step following spike (S)
del_G_syn_E_E = 9.5*10^(-9); %synaptic conductance step following spike (S)
del_G_syn_I_I = 1.4*del_G_syn_E_E; %synaptic conductance step following spike (S)
del_G_syn_E_I = 9.5*10^(-9); %synaptic conductance step following spike (S)
del_G_syn_I_E = 1.4*del_G_syn_E_E; %synaptic conductance step following spike (S)

%___________________________
del_G_sra = 330e-09; %spike rate adaptation conductance step following spike %ranges from 1-200 *10^(-9) (S)
tau_sra = (-2.4*10^5)*del_G_sra + 0.11; %spike rate adaptation time constant (s)
%If want to have STDP, change connectivity_gain to > 0.
connectivity_gain = 0; %0.005; %amount to increase or decrease connectivity by with each spike (more at the range of 0.002-0.005)

%How spikes are initiated:
%'cluster' sets a cluster to threshold;
%'current' means spikes depend on an input current of I_in; 
%'neuron' sets a random selection of 2% of neurons to threshold
type = 'current'; %'cluster'; %'neuron';

% %Rhythmic current input: (uncomment if desired)
% I_coeff = 0; %5.1*10^(-10); %set to 0 for no input current
% I_Hz = 1; %frequency of input - 1 Hz or 60 bpm, 4-8 Hz for theta rhythm
% I_in = I_coeff*(0.5 + 0.5*sin(I_Hz*2*pi*x_in)); %Approximately theta wave input current
%Noisy input: (uncomment if desired)
% I_coeff = 2.7; %2.7; %set to 0 for no input current
% I_scale = 1*10^(-9); %sets the scale of the current input
% Input conductance
G_coeff = -38; %-40;
G_scale = 1*10^(-9);

%Calculate connection probabilites
conn_prob = 0.08; %set a total desired connection probability
p_E = 0.75; %probability of an excitatory neuron

%Event Statistics
event_cutoff = 0.05; %0.25; %fraction of neurons that have to be involved to constitute a successful event

%Save parameters to a structure and to computer
w = whos;
parameters = struct;
for a = 1:length(w) 
    parameters.(w(a).name) = eval(w(a).name); 
end
clear w a
save(strcat(save_path,'/parameters.mat'),'parameters')

%% Create networks and test spike progressions
%Runs through a series of different random number generator seeds to change
%the network connectivity and setup, and then automatically outputs
%sequence data to a folder

%If uploading a parameter file, uncomment the next line
%load(strcat(save_path,'/parameters.mat'))

%____________________________________
%___Calculate Dependent Parameters___
%____________________________________
cluster_n = min(parameters.n*2/parameters.clusters,parameters.n); %number of neurons in a cluster (for small n round(n/3), for large n round(n/5)) 
parameters.('cluster_n') = cluster_n;

%Interaction constants
t_steps = parameters.t_max/parameters.dt; %number of timesteps in simulation
E_syn_E = parameters.V_syn_E*ones(parameters.n,1); %vector of the synaptic reversal potential for excitatory connections
E_syn_I = parameters.V_syn_I*ones(parameters.n,1); %vector of the synaptic reversal potential for inhibitory connections
%save for easy calculations
parameters.('t_steps') = t_steps;
parameters.('E_syn_E') = E_syn_E;
parameters.('E_syn_I') = E_syn_I;

%How many tests of different initializations to run
if strcmp(parameters.type,'cluster')
    test_val_max = parameters.clusters; %Ensures every cluster is initialized
else
    test_val_max = 10; %This value can be modified as you see fit
end
%save for easy calculations
parameters.('test_val_max') = test_val_max;

%Adding an input conductance to all cells (one of the options must be uncommented)
x_in = [0:parameters.dt:parameters.t_max];
% %Rhythmic current input: (uncomment if desired)
%Noisy input conductance: (uncomment if desired)
G_in = parameters.G_coeff*randn(parameters.n,parameters.t_steps+1)*parameters.G_scale;
%save for easy calculations
parameters.('x_in') = x_in;
parameters.('G_in') = G_in;

%Calculate connection probabilites
npairs = parameters.n*(parameters.n-1); %total number of possible neuron connections
nclusterpairs = parameters.cluster_n*(parameters.cluster_n - 1)*parameters.clusters; %total number of possible intra-cluster connections
cluster_prob = min(parameters.conn_prob*npairs/nclusterpairs,1); %0.2041; %intra-cluster connection probability
p_I = 1 - parameters.p_E; %probability of an inhibitory neuron
n_I = round(p_I*parameters.n); %number of inhibitory neurons
%save for easy calculations
parameters.('npairs') = npairs;
parameters.('nclusterpairs') = nclusterpairs;
parameters.('cluster_prob') = cluster_prob;
parameters.('p_I') = p_I;
parameters.('n_I') = n_I;

save(strcat(save_path,'/parameters.mat'),'parameters')

%____________________________________________
%___Run Simulations for Different Networks___
%____________________________________________

for i = 1:1%10 %how many different network structures to test
    rng(i) %set random number generator for network structure
    %CREATE NETWORK SAVE PATH
    net_save_path = strcat(save_path,'/network_',string(i));
    if ~isfolder(net_save_path)
        mkdir(net_save_path);
    end
    
    %SET UP NETWORK
    [cluster_mat, conns] = create_clusters(parameters, i, 1);
    conns_copy = conns; %just a copy of the connections to maintain for reset runs if there's "plasticity"
    %Randomize excitatory and inhibitory connection strengths based on selected
    %probability.
    all_indices = [1:parameters.n];
    I_indices = datasample(all_indices,parameters.n_I,'Replace',false); %indices of inhibitory neurons
    E_indices = find(~ismember(all_indices,I_indices)); %indices of excitatory neurons
    n_EE = sum(conns(E_indices,E_indices),'all'); %number of E-E connections
    n_EI = sum(conns(E_indices,I_indices),'all'); %number of E-I connections
    n_II = sum(conns(I_indices,I_indices),'all'); %number of I-I connections
    n_IE = sum(conns(I_indices,E_indices),'all'); %number of I-E connections
    clear all_indices
    
    %SAVE NETWORK STRUCTURE
    network = struct;
    network(1).cluster_mat = cluster_mat;
    network(1).conns = conns;
    network(1).I_indices = I_indices;
    network(1).E_indices = E_indices;
    save(strcat(net_save_path,'/network.mat'),'network')
    
    clear cluster_mat conns I_indices E_indices
    
    %RUN MODEL AND CALCULATE
    %Run through every cluster initialization and store relevant data and
    %calculations
    network_var = struct;
    network_spike_sequences = struct;
    network_cluster_sequences = struct; 
    
    for j = 1:parameters.test_val_max
        seed = j;
        
        %Create Storage Variables
        I_syn = zeros(parameters.n,parameters.t_steps+1); %synaptic current emitted by each neuron at each timestep (A)
        %synaptic conductances for each neuron at each timestep
        G_syn_I = zeros(parameters.n,parameters.t_steps+1); %conductance for presynaptic inhibitory (S)
        G_syn_E = zeros(parameters.n,parameters.t_steps+1); %conductance for presynaptic excitatory (S)
        V_m = zeros(parameters.n,parameters.t_steps+1); %membrane potential for each neuron at each timestep
        V_m(:,1) = parameters.V_reset + randn([parameters.n,1])*(10^(-3)); %set all neurons to baseline reset membrane potential with added noise
        G_sra = zeros(parameters.n,parameters.t_steps+1); %refractory conductance for each neuron at each timestep
        
        %Run model
        [V_m, G_sra, G_syn_I, G_syn_E, I_syn] = lif_sra_calculator_postrotation(...
            parameters, seed, network, I_syn, G_syn_I, G_syn_E, V_m, G_sra);
        network_var(j).V_m = V_m;
        network_var(j).G_sra = G_sra;
        network_var(j).G_syn_I = G_syn_I;
        network_var(j).G_syn_E = G_syn_E;
        network_var(j).I_syn = I_syn;
        
        %Find spike profile
        spikes_V_m = V_m >= parameters.V_th;
        [spikes_x,spikes_t] = find(spikes_V_m);
        max_time = max(spikes_t);
        spiking_neurons = unique(spikes_x, 'stable');
        
        %Calculate LFP
        