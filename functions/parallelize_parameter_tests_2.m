function avg_mat = parallelize_parameter_tests_2(parameters,num_nets,...
    num_inits, parameter_vec, test_n, ind)
    %_________
    %ABOUT: This function runs through a series of commands to test the
    %outputs of a particular parameter set in comparison to a strict set of
    %criteria. This function is to be used in conjunction with
    %network_tests2.m where a parallelized for loop calls this function,
    %parallelize_networks.m where network initializations are parallelized,
    %and parallelize_network_tests.m where network tests are parallelized.
    %
    %INPUTS:
    %   parameters = struct containing the following
    %       n = Number of neurons in the network
    %       V_reset = The reset membrane potential (V)
    %       V_m_noise = Magnitude of membrane potential simulation noise (V)
    %       V_th = The threshold membrane potential (V)
    %       
    %       del_G_sra = spike rate adaptation conductance step following spike 
    %               ranges from 1-200 *10^(-9) (S)
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
    %       E_syn_E = An [n x 1] vector of the synaptic reversal potential for
    %               excitatory connections (V)
    %       E_syn_I = An [n x 1] vector of the synaptic reversal potential for
    %               inhibitory connections (V)
    %       E_K = Potassium reversal potential (V)
    %       E_L = Leak reversal potential (V)
    %       G_L = Leak conductance (S) - 10-30 nS range
    %       C_m = Total membrane capacitance (F)
    %       dt = Timestep (s)
    %       tau_sra = Spike rate adaptation time constant (s)
    %       tau_syn_E = AMPA/NMDA synaptic decay time constant (s)
    %       tau_syn_I = GABA synaptic decay time constant (s)
    %       tau_stdp = STDP decay time constant (s)
    %       connectivity_gain = Amount to increase or decrease connectivity by 
    %               with each spike (more at the range of 1.002-1.005) -
    %               keep at 1 to ensure no connectivity change
    %       I_in = An [n x t_steps + 1] matrix with input current values
    %       t_steps = The number of timesteps in the simulation
    %       type = Determines how spiking is initiated. Either:
    %           1. type = 'cluster', which sets a cluster of neurons to
    %           threshold at step 1
    %           2. type = 'neuron', which sets a fraction of neurons to
    %           threshold at step 1, where the fraction is set on line 
    %           3. type = 'current', which means spiking depends entirely on
    %           I_in and noise for initialization
    %   num_nets = number of network structures to test per parameter set 
    %   num_inits = number of network initializations to test per network
    %       structure
    %   parameter_vec = 
    %OUTPUTS:
    %   avg_mat = A vector representing average values from all network
    %   initializations and the firing initializations for each network
    %   structure. Namely, the first dimension length is the number of
    %   networks, the second dimension length is the number of
    %   initializations, and the third dimension length is the 3 parameters
    %   from each combination:
    %       1. number of spiking neurons
    %       2. average firing rate
    %       3. average event length
    %_________
    %Get index values for parameter combination
    [num_params, ~] = size(parameter_vec);
    indices = 1:test_n^num_params;
    size_vec = test_n * ones(1,num_params);
    [ind_1,ind_2] = ind2sub(size_vec,indices);
    %Pull parameter combination
    parameters.del_G_sra = parameter_vec(1,ind_1(ind));
    parameters.tau_sra = (-2.4*10^5)*parameters.del_G_sra + 0.11;
    parameters.G_coeff = parameter_vec(2,ind_2(ind));
    G_in = parameters.G_coeff*randn(parameters.n,parameters.t_steps+1)*parameters.G_scale;
    parameters.G_in = G_in;
    %Run network initialization code
    resp_mat = zeros(num_nets, 3);
    parfor i = 1:num_nets, resp_mat(i,:) = parallelize_networks_2(parameters,i, num_inits); end
    avg_mat = mean(resp_mat,1);
    disp(strcat('Index #',string(ind),'= complete'))
end