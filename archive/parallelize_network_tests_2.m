function [outputVec, allResults] = parallelize_network_tests_2(parameters, network, j)
    %_________
    %ABOUT: This function runs through a series of commands to test the
    %outputs of a particular parameter set in comparison to a strict set of
    %criteria. This function is to be used in conjunction with
    %network_tests2.m where a parallelized for loop calls this function.
    %
    %INPUTS:
    %   parameters = a structure that contains the following (only
    %   relevant listed below):
    %       n = Number of neurons in the network
    %       clusters = Number of clusters of neurons in network
    %       t_max = maximum time of simulation (s)
    %       dt = timestep of simulation (s)
    %       tau_syn_E = AMPA synaptic decay time constant (s) [Ignoring NMDA as slow and weak]
    %       tau_syn_I = GABA synaptic decay time constant (s)
    %       tau_stdp = STDP decay time constant (s)
    %       E_K = Potassium reversal potential (V)
    %       E_L = Leak reversal potential (V)
    %       G_L = Leak conductance (S) - 10-30 nS range
    %       C_m = Total membrane capacitance (F)
    %       V_m_noise = Magnitude of membrane potential simulation noise (V)
    %       V_th = The threshold membrane potential (V)
    %       V_reset = The reset membrane potential (V)
    %       V_syn_E = Excitatory synaptic reversal potential (V)
    %       V_syn_I = Inhibitory synaptic reversal potential (V)
    %       del_G_syn_E_E = Synaptic conductance step for 
    %               excitatory-excitatory neuron connections following
    %               spike (S)
    %       del_G_syn_E_I = Synaptic conductance step for 
    %               excitatory-inhibitory neuron connections following
    %               spike (S)
    %       del_G_syn_I_I = Synaptic conductance step for 
    %               inhibitory-inhibitory neuron connections following
    %               spike (S)
    %       del_G_syn_I_E = Synaptic conductance step for 
    %               inhibitory-excitatory neuron connections following
    %               spike (S)
    %       del_G_sra = spike rate adaptation conductance step following spike 
    %               ranges from 1-200 *10^(-9) (S)
    %       tau_sra = Spike rate adaptation time constant (s)
    %       connectivity_gain = Amount to increase or decrease connectivity by 
    %               with each spike (more at the range of 1.002-1.005) -
    %               keep at 1 to ensure no connectivity change
    %       G_coeff = input conductance coefficient (setting strength) (S)
    %       G_scale = input conductance scale (ex. nano = 1*10^(-9)) (S)
    %       t_steps = The number of timesteps in the simulation
    %       syn_E = An [n x 1] vector of the synaptic reversal potential for
    %               excitatory connections (V)
    %       syn_I = An [n x 1] vector of the synaptic reversal potential for
    %               inhibitory connections (V)     
    %   network = a structure that contains the following:
    %       cluster_mat = A binary [clusters x n] matrix of which neurons are
    %               in which cluster
    %       conns = An [n x n] matrix of which neurons are connected to each
    %               other, with values greater than 1 implying stronger
    %               connectivity
    %       I_indices = Vector of indices of inhibitory neurons
    %       E_indices = Vector of indices of excitatory neurons
    %   j = random number generator seed
    %OUTPUTS:
    %   vec = A vector of parameters for a run:
    %       1. number of spiking neurons
    %       2. average firing rate
    %       3. average event length
    %_________
    
    seed = j;
    
    %Create input conductance variable
    if parameters.usePoisson
        G_in = zeros(parameters.n, parameters.t_steps+1);
        for k = 2:(parameters.t_steps+1)
            G_in(:,k) = G_in(:,k-1)*exp(-parameters.dt/parameters.tau_syn_E);
            G_in(:,k) = G_in(:,k) + parameters.W_gin * [rand(parameters.n, 1) < (parameters.dt*parameters.rG)];
        end
    else
        G_in = (parameters.G_std*randn(parameters.n,parameters.t_steps+1))+parameters.G_mean;
        G_in(G_in<0) = 0;
    end
    parameters.('G_in') = G_in;
        
    %Create Storage Variables
    V_m = zeros(parameters.n,parameters.t_steps+1); %membrane potential for each neuron at each timestep
    V_m(:,1) = parameters.V_reset + randn([parameters.n,1])*parameters.V_m_noise; %set all neurons to baseline reset membrane potential with added noise
    %Run model
    [V_m, ~, ~, ~, ~] = randnet_calculator(parameters, seed, network, V_m);
    clear I_syn G_syn_I G_syn_E
    
    %Find spike profile
    spikes_V_m = V_m >= parameters.V_th;
    [spikes_x,spikes_t] = find(spikes_V_m);
    spiking_neurons = unique(spikes_x, 'stable');
    
    % detect events and compute outputs
    network_spike_sequences = struct; 
    network_cluster_sequences = struct;
    [network_spike_sequences, network_cluster_sequences, outputVec] = detect_events(parameters, network, V_m , j, network_spike_sequences, network_cluster_sequences);
    
    % Overall simulation statistics
    allResults.ithInit = j;
    try
        allResults.numEvents = numel(network_spike_sequences(j).event_lengths); % number of detected events
    catch
        allResults.numEvents = [];
    end
    allResults.fracFire =  mean(sum(spikes_V_m, 2)>0); % Fraction of cells that fire at all during simulation
    allResults.meanRate = mean(sum(spikes_V_m, 2)/parameters.t_max); % mean over cells' average firing rate
    allResults.stdRate = std(sum(spikes_V_m, 2)/parameters.t_max); % STD over cells' average firing rate
    
    % Stats for each detected event
    try
        allResults.eventLength = network_spike_sequences(j).event_lengths; % duration in seconds of all detected events
    catch
        allResults.eventLength = [];
    end
    try
        allResults.eventParticipation = structfun( @mean , network_spike_sequences(j).nonspiking_neurons )'; % fraction of cells that fired in each event
    catch
        allResults.eventParticipation = [];
    end
    
end