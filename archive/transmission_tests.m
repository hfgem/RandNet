%Test how activity propagates through the network

%% Test how a single neuron spike affects connected neurons

%Load network structure and parameters
% load_path = uigetdir('/Users/hannahgermaine/Documents/PhD/','Select Network Load Folder'); %Have user input where they're pulling parameters from
% load(strcat(load_path,'/network.mat'))
% slashes = find(load_path == '/');
% param_path = load_path(1:slashes(end));
% load(strcat(param_path,'/parameters.mat'))

%Change simulation time and initialization type
parameters_copy.t_max = 0.1;
parameters_copy.type = 'neuron';

%Select how many simulation steps following a spike to observe
view_steps = 10;

%________________________________________
%______Calculate dependent variables______
%_________________________________________

%Interaction constants
t_steps = parameters_copy.t_max/parameters_copy.dt; %number of timesteps in simulation
E_syn_E = parameters_copy.V_syn_E*ones(parameters_copy.n,1); %vector of the synaptic reversal potential for excitatory connections
E_syn_I = parameters_copy.V_syn_I*ones(parameters_copy.n,1); %vector of the synaptic reversal potential for inhibitory connections
IES = ceil(parameters_copy.IEI/parameters_copy.dt); %inter-event-steps = the number of steps to elapse between spikes
%save for easy calculations
parameters_copy.('t_steps') = t_steps;
parameters_copy.('E_syn_E') = E_syn_E;
parameters_copy.('E_syn_I') = E_syn_I;
parameters_copy.('IES') = IES;

%How many tests of different initializations to run
if strcmp(parameters_copy.type,'cluster')
    test_val_max = parameters_copy.clusters; %Ensures every cluster is initialized
else
    test_val_max = 10; %This value can be modified as you see fit
end
%save for easy calculations
parameters_copy.('test_val_max') = test_val_max;

%Adding an input current to all cells (one of the options must be uncommented)
x_in = [0:parameters_copy.dt:parameters_copy.t_max];
% %Rhythmic current input: (uncomment if desired)
% I_coeff = 0; %5.1*10^(-10); %set to 0 for no input current
% I_Hz = 1; %frequency of input - 1 Hz or 60 bpm, 4-8 Hz for theta rhythm
% I_in = I_coeff*(0.5 + 0.5*sin(I_Hz*2*pi*x_in)); %Approximately theta wave input current
%Noisy input: (uncomment if desired)
I_in = parameters_copy.I_coeff*randn(parameters_copy.n,parameters_copy.t_steps+1)*parameters_copy.I_scale; %Generally use -0.5-0.5 nA stimulus
%save for easy calculations
parameters_copy.('x_in') = x_in;
parameters_copy.('I_in') = I_in;

%Calculate connection probabilites
npairs = parameters_copy.n*(parameters_copy.n-1); %total number of possible neuron connections
nclusterpairs = parameters_copy.cluster_n*(parameters_copy.cluster_n - 1)*parameters_copy.clusters; %total number of possible intra-cluster connections
cluster_prob = min(parameters_copy.conn_prob*npairs/nclusterpairs,1); %0.2041; %intra-cluster connection probability
p_I = 1 - parameters_copy.p_E; %probability of an inhibitory neuron
n_I = round(p_I*parameters_copy.n); %number of inhibitory neurons
%save for easy calculations
parameters_copy.('npairs') = npairs;
parameters_copy.('nclusterpairs') = nclusterpairs;
parameters_copy.('cluster_prob') = cluster_prob;
parameters_copy.('p_I') = p_I;
parameters_copy.('n_I') = n_I;

%% Run simulation and look at post-spike deflections
for i = 1:10
    %Load network parameters
    seed = i;
    
    %CREATE SAVE PATH
    net_save_path = strcat(load_path,'/deflection_tests');
    if ~isfolder(net_save_path)
        mkdir(net_save_path);
    end
    
    %Create Storage Variables
    I_syn = zeros(parameters_copy.n,parameters_copy.t_steps+1); %synaptic current emitted by each neuron at each timestep (A)
    %synaptic conductances for each neuron at each timestep
    G_syn_I = zeros(parameters_copy.n,parameters_copy.t_steps+1); %conductance for presynaptic inhibitory (S)
    G_syn_E = zeros(parameters_copy.n,parameters_copy.t_steps+1); %conductance for presynaptic excitatory (S)
    V_m = zeros(parameters_copy.n,parameters_copy.t_steps+1); %membrane potential for each neuron at each timestep
    V_m(:,1) = parameters_copy.V_reset + randn([parameters_copy.n,1])*(10^(-3)); %set all neurons to baseline reset membrane potential with added noise
    G_sra = zeros(parameters_copy.n,parameters_copy.t_steps+1); %refractory conductance for each neuron at each timestep
    
    %Now run simulation
    [V_m, G_sra, G_syn_I, G_syn_E, I_syn] = lif_sra_calculator_postrotation(...
            parameters_copy, seed, network, I_syn, G_syn_I, G_syn_E, V_m, G_sra);
    
    %Pull out post-synaptic neuron parameters
    first_spiker_index = find(V_m(:,1) >= parameters_copy.V_th);
    first_spiker_index = first_spiker_index(1);
    connected_neurons = find(network.conns(first_spiker_index,:));
    V_m_connected = V_m(connected_neurons,1:view_steps); %membrane potential
%     G_syn_I_connected = G_syn_I(connected_neurons,1:view_steps); %inhibitory conductance
%     G_syn_E_connected = G_syn_E(connected_neurons,1:view_steps); %excitatory conductance
    
    %Calculate average and variance of post-synaptic parameters
    V_m_avg = mean(V_m_connected,1)*10^3; %converted to mV scale
    V_m_var = std(V_m_connected,1)*10^3; %converted to mV scale
%     G_syn_I_avg = mean(G_syn_I_connected,1);
%     G_syn_I_var = std(G_syn_I_connected,1);
%     G_syn_E_avg = mean(G_syn_E_connected,1);
%     G_syn_E_var = std(G_syn_E_connected,1);
    
    %Plot
    f = figure;
    %Membrane potential
    inBetween = [V_m_avg + V_m_var, fliplr(V_m_avg - V_m_var)];
    fill([0:view_steps-1, fliplr(0:view_steps-1)], inBetween, 'b','FaceAlpha',0.5);
    hold on
    plot(0:view_steps-1,V_m_avg)
    for j = 1:length(connected_neurons)
        scatter([0:view_steps - 1],V_m_connected(j,:)*10^(3),'filled')
    end
    xlabel('Timesteps Following Spike')
    ylabel('Membrane Potential (mV)')
    title(strcat('V_m Deflections'))
    %Plot Title
    sgtitle('Average and Variance Deflections From Single Neuron Spike')
    
    %Save Figure
    savefig(f,strcat(net_save_path,'/deflections_from_neuron_',string(first_spiker_index),'.fig'))
    saveas(f,strcat(net_save_path,'/deflections_from_neuron_',string(first_spiker_index),'.jpg'))
    
    close(f)
end    

%% 