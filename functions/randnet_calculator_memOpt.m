function [spikeMat, conns] = randnet_calculator_memOpt(parameters, seed, network, V_m)
    %_________
    % ABOUT: This function produces identical results as
    % randnet_calculator, but the only variables it saves is V_m.
    % Use randnet_calculator for simulations where you want to analyze
    % other variables. 
    % Use randnet_calculator for simulations where you only need V_m and
    % wanted to reduce memory usage. 
    % 
    % The script test_randnet_calculator tests that randnet_calcualtor and
    % randnet_calculator_memOpt produce the same V_m and conns
    %
    % Note: another option is to have a single simulation function with a 
    % flag inside the loop that can save the variable vectors at each time
    % step. But this if statement slows the simulation by ~3%
    
    
    %Set the random number generator seed
    rng(seed)
    

    %Create Storage Variables (single time-step, updated each increment)
    G_sra = zeros(parameters.n,1); %refractory conductance for each neuron at each timestep (S)
    G_syn_I_E = zeros(parameters.n,1); %conductance for pre-inhib to post-excit (S)
    G_syn_E_E = zeros(parameters.n,1); %conductance for pre-excit to post-excit (S)
    G_syn_I_I = zeros(parameters.n,1); %conductance for pre-inhib to post-inhib (S)
    G_syn_E_I = zeros(parameters.n,1); %conductance for pre-excit to post-inhib (S)

    %Copy connectivity matrix in case of stdp changes
    conns = network.conns; %separately update a connectivity matrix
    
    %Binary indices of excitatory and inhibitory neurons
    E_bin = zeros(parameters.n,1);
    E_bin(network.E_indices) = 1;
    I_bin = zeros(parameters.n,1);
    I_bin(network.I_indices) = 1;
    
    %Variables for STDP
    t_spike = zeros(parameters.n,1); %vector to store the time of each neuron's last spike, for use in STDP
    t_stdp = round(parameters.tau_stdp/parameters.dt);
    
    spikeMat = false(parameters.n, parameters.t_steps+1); 
    
    %Run through each timestep and calculate
    for t = 1:parameters.t_steps
        %check for spiking neurons and postsynaptic and separate into E and I
        spikers = find(V_m >= parameters.V_th);
        spikeMat(spikers,t) = spikers;
        
        t_spike(spikers) = t;
        spikers_I = spikers(ismember(spikers,network.I_indices)); %indices of inhibitory spiking presynaptic neurons
        spikers_E = spikers(ismember(spikers,network.E_indices)); %indices of excitatory spiking presynaptic neurons
        %______________________________________
        %Adjust parameters dependent on spiking
        G_sra(spikers) = G_sra(spikers) + parameters.del_G_sra; %set SRA conductance values
        %Synaptic conductance is stepped for postsynaptic neurons
        %   dependent on the number of presynaptic connections, and the
        %   current will depend on the presynaptic neuron type (E_syn_I and E_syn_E)
        incoming_conn_E = sum(conns(spikers_E,:),1)'; %post-synaptic neuron E input counts
        incoming_conn_I = sum(conns(spikers_I,:),1)'; %post-synaptic neuron I input counts
        G_syn_I_E = G_syn_I_E + parameters.del_G_syn_I_E*incoming_conn_I.*E_bin;
        G_syn_E_E = G_syn_E_E + parameters.del_G_syn_E_E*incoming_conn_E.*E_bin;
        G_syn_I_I = G_syn_I_I + parameters.del_G_syn_I_I*incoming_conn_I.*I_bin;
        G_syn_E_I = G_syn_E_I + parameters.del_G_syn_E_I*incoming_conn_E.*I_bin;
        %______________________________________
        %Calculate membrane potential using integration method
        V_ss = ( parameters.G_in(:,t).*parameters.syn_E + G_syn_E_E.*parameters.syn_E + G_syn_E_I.*parameters.syn_E + G_syn_I_I.*parameters.syn_I + G_syn_I_E.*parameters.syn_I + parameters.G_L*parameters.E_L + G_sra*parameters.E_K )./(parameters.G_L + G_sra + G_syn_E_E + G_syn_E_I + G_syn_I_I + G_syn_I_E + parameters.G_in(:,t));
        taueff = parameters.C_m./(parameters.G_L + G_sra + G_syn_E_E + G_syn_E_I + G_syn_I_I + G_syn_I_E + parameters.G_in(:,t));
        V_m = V_ss + (V_m - V_ss).*exp(-parameters.dt ./taueff) + randn([parameters.n,1])*parameters.V_m_noise*sqrt(parameters.dt); %the randn portion can be removed if you'd prefer no noise
        V_m(spikers) = parameters.V_reset; %update those that just spiked to reset
        %______________________________________
        %Update next step conductances
        G_sra = G_sra*exp(-parameters.dt/parameters.tau_sra); %Spike rate adaptation conductance
        %Synaptic conductance updated for each postsynaptic neuron by
        %incoming connection type
        G_syn_E_E = G_syn_E_E.*exp(-parameters.dt/parameters.tau_syn_E); %excitatory conductance update
        G_syn_I_E = G_syn_I_E.*exp(-parameters.dt/parameters.tau_syn_I); %excitatory conductance update
        G_syn_I_I = G_syn_I_I.*exp(-parameters.dt/parameters.tau_syn_I); %inhibitory conductance update
        G_syn_E_I = G_syn_E_I.*exp(-parameters.dt/parameters.tau_syn_E); %inhibitory conductance update
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