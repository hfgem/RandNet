function pass = parallelize_network_tests(parameters,network,j)
    %_________
    %ABOUT: This function runs through a series of commands to test the
    %outputs of a particular parameter set in comparison to a strict set of
    %criteria. This function is to be used in conjunction with
    %network_tests2.m where a parallelized for loop calls this function.
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
    %   pass = A binary value of whether the parameters resulted in output
    %       that passed the criteria.
    %_________
    
    pass = 0;
    seed = j;
    %Create Storage Variables
    I_syn = zeros(parameters.n,parameters.t_steps+1); %synaptic current emitted by each neuron at each timestep (A)
    %synaptic conductances for each neuron at each timestep
    G_syn_I = zeros(parameters.n,parameters.t_steps+1); %conductance for presynaptic inhibitory (S)
    G_syn_E = zeros(parameters.n,parameters.t_steps+1); %conductance for presynaptic excitatory (S)
    V_m = zeros(parameters.n,parameters.t_steps+1); %membrane potential for each neuron at each timestep
    V_m(:,1) = parameters.V_reset + randn([parameters.n,1])*parameters.V_m_noise; %set all neurons to baseline reset membrane potential with added noise
    G_sra = zeros(parameters.n,parameters.t_steps+1); %refractory conductance for each neuron at each timestep
    %Run model
    [V_m, ~, ~, ~, ~] = lif_sra_calculator_postrotation(...
        parameters, seed, network, I_syn, G_syn_I, G_syn_E, V_m, G_sra);
    clear I_syn G_syn_I G_syn_E
    
    %Find spike profile
    spikes_V_m = V_m >= parameters.V_th;
    [spikes_x,spikes_t] = find(spikes_V_m);
    spiking_neurons = unique(spikes_x, 'stable');
    %First check if the number of spiking neurons fits criteria
    if length(spiking_neurons) >= parameters.event_cutoff*parameters.n
        %Find maximum firing rate + average maximum firing rates of neurons
        all_fr = sum(spikes_V_m(spiking_neurons,:),2)/parameters.t_max; %only look at neurons that fired at least 1x
        avg_fr = mean(all_fr);
        %Test of success (passing_trials = passing_trials + 1;)
        %According to Hirase et al. 2001, neurons fire
        %around 1 Hz on average and have a max around 1.5 Hz
        if and(avg_fr>= 0.02, avg_fr <= 1.5)
            %Find firing event times
            events = [];
            event_lengths = [];
            last_start = spikes_t(1);
            last_time = spikes_t(1);
            spike_count = 1;
            for t_i = 2:length(spikes_t)
                s_i = spikes_t(t_i);
                if s_i - last_time <= parameters.IES
                    last_time = s_i;
                    spike_count = spike_count + 1;
                else
                    if (last_start ~= last_time) && (spike_count > parameters.event_cutoff*parameters.n) %weed out events w/ too few spikes
                        events = [events; [last_start, last_time]]; %add the last range of spikes to the events vector
                        event_lengths = [event_lengths, (last_time - last_start)*parameters.dt]; %#ok<*AGROW>
                    end
                    last_start = s_i;
                    last_time = s_i;
                    spike_count = 1;
                end
            end
            clear t_i s_i
            if (last_start ~= last_time) && (spike_count > parameters.event_cutoff*parameters.n) %weed out events w/ too few spikes
                events = [events; [last_start, last_time]]; %add the last interval
                event_lengths = [event_lengths, (last_time - last_start)*parameters.dt]; %#ok<*SAGROW>
            end
            avg_event_length = mean(event_lengths);
            %Ylinen et al 1995 - average SWR event lasts between 50 ms and
            %120 ms - we extended the boundaries by 30 ms
            if and(avg_event_length >= 0.02, avg_event_length <= 0.15)
                pass = 1;
%                 %Check membrane potentials
%                 %check first event for membrane
%                 %potential (good enough!)
%                 event_V_m = V_m(spiking_neurons,events(1,1):events(1,2));
%                 std_V_m = std(event_V_m,1,'all');
%                 %Ylinen et al. 1995 1-5 mV V_m oscillation during SWR
%                 if and(std_V_m >= 1*10^(-3), std_V_m <= 5*10^(-3))
%                     pass = 1; %Success! Passed all strict criteria!
%                 end %membrane potential if loop  
            end %SWR length if loop
        end %firing rate if loop
    end %number of neurons if loop
end