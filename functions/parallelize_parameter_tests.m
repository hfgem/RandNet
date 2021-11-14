function success = parallelize_parameter_tests(parameters,num_nets,...
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
    %OUTPUTS:
    %   success = A value in the range of [0,1] describing the fraction of
    %       successful initializations (passing criteria).
    %_________
    %Get index values for parameter combination
    indices = 1:test_n^4;
    [tau_ind,del_G_sra_ind,del_G_syn_E_ind,del_G_syn_I_ind] = ind2sub([test_n,test_n,test_n,test_n],indices);
    %Pull parameter combination
    parameters.tau_sra = parameter_vec(1,tau_ind(ind));
    parameters.del_G_sra = parameter_vec(2,del_G_sra_ind(ind));
    parameters.del_G_syn_E_E = parameter_vec(3,del_G_syn_E_ind(ind));
    parameters.del_G_syn_E_I = parameter_vec(3,del_G_syn_E_ind(ind));
    parameters.del_G_syn_I_E = parameter_vec(4,del_G_syn_I_ind(ind));
    parameters.del_G_syn_I_I = parameter_vec(4,del_G_syn_I_ind(ind));
    %Run network initialization code
    passing_trials = zeros(1,num_nets);
    parfor i = 1:num_nets, passing_trials(i) = parallelize_networks(parameters,i, num_inits); end
    passing_trials = sum(passing_trials,'all');
    success = passing_trials/(num_nets*num_inits); %store fraction of successful trials
    disp([ind, success])
end