function [V_m, G_sra, G_syn_I, G_syn_E, conns] = lif_sra_calculator_postrotation(...
    parameters, seed, network, G_syn_I, G_syn_E, V_m, G_sra)
    %_________
    %ABOUT: This function uses the leaky integrate-and-fire model of 
    %neuronal firing to calculate the trajectory of membrane potentials,
    %currents, etc... that take place in a particular network with a 
    %particular set of parameters and initialization.
    %
    %INPUTS:
    %   parameters = a structure that contains the following:
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
    %   seed = A random number generator seed which:
    %       1. when type = 'cluster' sets which cluster is to be used for
    %           the initialization of spiking
    %       2. when type = 'neuron' sets the random number generator seed
    %           which affects the random selection of neurons used in the
    %           initialization of spiking as well as the random noise added
    %           to the membrane potential
    %       3. when type = 'current' sets the random number generator seed 
    %           which affects only the random noise added to the membrane 
    %           potential 
    %   network = a structure that contains the following:
    %       cluster_mat = A binary [clusters x n] matrix of which neurons are
    %               in which cluster
    %       conns = An [n x n] matrix of which neurons are connected to each
    %               other, with values greater than 1 implying stronger
    %               connectivity
    %       I_indices = Vector of indices of inhibitory neurons
    %       E_indices = Vector of indices of excitatory neurons
    %   G_syn_I = An [n x t_steps+1] matrix of conductance for presynaptic 
    %               inhibitory (S)
    %   G_syn_E = An [n x t_steps+1] matrix of conductance for presynaptic 
    %               excitatory (S)
    %   V_m = An [n x t_steps+1] matrix of membrane potential for each 
    %               neuron at each timestep
    %   G_sra = An [n x t_steps+1] matrix with refractory conductance for
    %               each neuron at each timestep (S)
    %
    %OUTPUTS:
    %   V_m = An [n x t_steps+1] matrix of membrane potential for each 
    %               neuron at each timestep
    %   G_sra = An [n x t_steps+1] matrix with refractory conductance for
    %               each neuron at each timestep (S)
    %   G_syn_I = An [n x t_steps+1] matrix of conductance for presynaptic 
    %               inhibitory (S)
    %   G_syn_E = An [n x t_steps+1] matrix of conductance for presynaptic 
    %               excitatory (S)
    %   conns = The updated conns matrix (if connectivity_gain != 1)
    %
    %ASSUMPTIONS:
    %   1. LIF model with spike rate adaptation and synaptic transmission
    %   2. SRA has decay
    %   3. Synaptic transmission has decay
    %   4. There can be an externally applied current input through 
    %       I_in
    %   5. Excitatory and inhibitory connections have different reversal
    %       potentials in the postsynaptic neuron represented in vectors
    %       E_syn_E and E_syn_I
    %   6. Excitatory and inhibitory currents have different synaptic
    %       conductance steps and decay rates
    %_________

    %Pick which cluster is initially spiking
    if strcmp(parameters.type,'cluster') %cluster set to threshold at timestep 1
        rng(1)
        neur_start = find(network.cluster_mat(seed,:)); %one cluster starts firing
        %Set the membrane potential of spiking neurons to threshold
        V_m(neur_start,1) = parameters.V_th; %#ok<FNDSB>
    elseif strcmp(parameters.type,'neuron')
        rng(seed)
        neur_start = rand(parameters.n,1) <= 0.05;
        %Set the membrane potential of spiking neurons to threshold
        V_m(neur_start,1) = parameters.V_th; 
    elseif strcmp(parameters.type,'current')
        rng(seed)
    end
    
    G_in = parameters.G_coeff*randn(parameters.n,parameters.t_steps+1)*parameters.G_scale;
    parameters.('G_in') = G_in; %save for easy calculations
    
    conns = network.conns; %separately update a connectivity matrix
    
    %TO DEPRECATE LATER: In the case that the parameters file being used 
    %does not contain V_m_noise:
    try
        parameters.V_m_noise;
    catch
        parameters.('V_m_noise') = 10^(-4);
    end
    
    %TO DEPRECATE LATER: In the case that the parameters file being used 
    %does not contain tau_stdp:
    try
        parameters.tau_stdp;
    catch
        parameters.('tau_stdp') = 5*10^(-3);
    end
    
    %Binary indices of excitatory and inhibitory neurons
    E_bin = zeros(parameters.n,1);
    E_bin(network.E_indices) = 1;
    I_bin = zeros(parameters.n,1);
    I_bin(network.I_indices) = 1;
    
    %Variables for STDP
    t_spike = zeros(parameters.n,1); %vector to store the time of each neuron's last spike, for use in STDP
    t_stdp = round(parameters.tau_stdp/parameters.dt);
    
    %Run through each timestep and calculate
    for t = 1:parameters.t_steps
        %check for spiking neurons and postsynaptic and separate into E and I
        spikers = find(V_m(:,t) >= parameters.V_th);
        t_spike(spikers) = t;
        spikers_I = spikers(ismember(spikers,network.I_indices)); %indices of inhibitory presynaptic neurons
        spikers_E = spikers(ismember(spikers,network.E_indices)); %indices of excitatory presynaptic neurons
        %______________________________________
        %Adjust parameters dependent on spiking
        G_sra(spikers,t) = G_sra(spikers,t) + parameters.del_G_sra; %set SRA conductance values
        %Synaptic conductance is stepped for postsynaptic neurons
        %   dependent on the number of presynaptic connections, and the
        %   current will depend on the presynaptic neuron type (E_syn_I and E_syn_E)
        incoming_conn_E = sum(conns(spikers_E,:),1)'; %post-synaptic neuron E input counts
        incoming_conn_I = sum(conns(spikers_I,:),1)'; %post-synaptic neuron I input counts
        G_syn_E(:,t) = G_syn_E(:,t) + parameters.del_G_syn_E_E*incoming_conn_E.*E_bin + parameters.del_G_syn_E_I*incoming_conn_E.*I_bin; %excitatory conductance
        G_syn_I(:,t) = G_syn_I(:,t) + parameters.del_G_syn_I_I*incoming_conn_I.*I_bin + parameters.del_G_syn_I_E*incoming_conn_I.*E_bin; %inhibitory conductance
        %______________________________________
        %Calculate membrane potential using integration method
        V_ss = ( parameters.G_in(:,t).*parameters.E_syn_E + G_syn_E(:,t).*parameters.E_syn_E + G_syn_I(:,t).*parameters.E_syn_I + parameters.G_L*parameters.E_L + G_sra(:,t)*parameters.E_K )./(parameters.G_L + G_sra(:,t) + G_syn_E(:,t) + G_syn_I(:,t) + parameters.G_in(:,t));
        taueff = parameters.C_m./(parameters.G_L + G_sra(:,t) + G_syn_E(:,t) + G_syn_I(:,t) + parameters.G_in(:,t));
        V_m(:,t+1) = V_ss + (V_m(:,t) - V_ss).*exp(-parameters.dt ./taueff) + randn([parameters.n,1])*parameters.V_m_noise; %the randn portion can be removed if you'd prefer no noise
        V_m(spikers,t+1) = parameters.V_reset; %update those that just spiked to reset
        %______________________________________
        %Update next step conductances
        G_sra(:,t+1) = G_sra(:,t)*exp(-parameters.dt/parameters.tau_sra); %Spike rate adaptation conductance
        %Synaptic conductance updated for each postsynaptic neuron by
        %incoming connection type
        G_syn_E(:,t+1) = G_syn_E(:,t).*exp(-parameters.dt/parameters.tau_syn_E); %excitatory conductance update
        G_syn_I(:,t+1) = G_syn_I(:,t).*exp(-parameters.dt/parameters.tau_syn_I); %inhibitory conductance update
        %______________________________________
        %Update connection strengths via STDP
        pre_syn_n = sum(conns(:,spikers),2) > 0; %pre-synaptic neurons
        post_syn_n = sum(conns(spikers,:),1) > 0; %post-synaptic neurons
        pre_syn_t = t_spike.*pre_syn_n; %spike times of pre-synaptic neurons
        post_syn_t = t_spike.*post_syn_n'; %spike times of post-synaptic neurons
        t_diff_pre = t - pre_syn_t; %time diff between pre-synaptic and current
        t_diff_post = t - post_syn_t; %time diff between post-synaptic and current
        del_conn_pre = parameters.connectivity_gain*exp(-t_diff_pre/t_stdp);
        del_conn_post = parameters.connectivity_gain*exp(-t_diff_post/t_stdp);
        conns(:,spikers) = conns(:,spikers) + del_conn_pre - del_conn_post; %enhance connections of those neurons that just fired
    end
end