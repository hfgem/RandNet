function avg_mat = parallelize_networks_2(parameters, i, num_inits, save_path)
    %_________
    %ABOUT: This function runs through a series of commands to test the
    %outputs of a particular parameter set in comparison to a strict set of
    %criteria. This function is to be used in conjunction with
    %network_tests2.m where a parallelized for loop calls this function.
    %This function then calls parallelize_network_tests.m to further
    %parallelize the code.
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
    %   i = random number generator seed
    %   num_inits = number of initializations to test
    %OUTPUTS:
    %   avg_mat = A vector of averaged resulting values from different
    %   initializations. Namely, the rows are different initializations,
    %   and the columns are:
    %       1. number of spiking neurons
    %       2. average firing rate
    %       3. average event length
    %_________

    rng(i,'twister') %set random number generator for network structure
    all_indices = [1:parameters.n];
    I_indices = datasample(all_indices,parameters.n_I,'Replace',false); %indices of inhibitory neurons
    E_indices = find(~ismember(all_indices,I_indices)); %indices of excitatory neurons
    [cluster_mat, conns] = create_clusters(parameters, i, 1);
    clear all_indices
    %Add in global inhibition, added to individual connections already
    %given. If global inhibition overrides any pre-set connections with
    %inhibitory neurons, reset values to global inhibition values.
    if parameters.only_global
        conns(I_indices,:) = parameters.I_strength*(rand([parameters.n*(1-parameters.p_E), parameters.n]) < parameters.p_I);
    else
        conns(I_indices,:) = conns(I_indices,:) + parameters.I_strength*(rand([parameters.n*(1-parameters.p_E), parameters.n]) < parameters.p_I);
    end
    %Save network structure
    network = struct;
    network(1).cluster_mat = cluster_mat;
    network(1).conns = conns;
    network(1).I_indices = I_indices;
    network(1).E_indices = E_indices;
    clear cluster_mat conns I_indices E_indices
    mat = zeros(num_inits,3);
    parfor j = 1:num_inits, mat(j,:) = parallelize_network_tests_2(parameters, ...
            network,j, save_path); end
    mat(isnan(mat)) = 0;
    avg_mat = sum(mat,1) ./ sum(mat > 0,1); %Only averaging those that did successfully produce data
end