function vec = parallelize_network_tests_2(parameters,network,j, save_path)
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
    
    vec = zeros(1,3);
    seed = j;
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
    
    %First value: the number of spiking neurons
    vec(1) = length(spiking_neurons);
    
    if vec(1) > 0
        %Find maximum firing rate + average maximum firing rates of neurons
        all_fr = sum(spikes_V_m(spiking_neurons,:),2)/parameters.t_max; %only look at neurons that fired at least 1x
        avg_fr = mean(all_fr);
        %Test of success (passing_trials = passing_trials + 1;)
        %According to Hirase et al. 2001, neurons fire
        %around 1 Hz on average and have a max around 1.5 Hz

        %Second value: the average firing rate
        vec(2) = avg_fr;
        
        if vec(2) > 0
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
            
            [num_events,~] = size(events);
            
            %Select results to visualize and save
            if avg_event_length > 0
                %if length(spiking_neurons) >= parameters.event_cutoff*parameters.n
                    %if and(avg_fr>= 0.02, avg_fr <= 1.5)
                        %if and(avg_event_length >= 0.02, avg_event_length <= 0.15)
                            %to_vis = rand;
                            %if to_vis <= 0.01
                                f = figure;
                                for e_i = 1:num_events
                                    subplot(1,num_events,e_i)
                                    imagesc(spikes_V_m(spiking_neurons,events(e_i,1):events(e_i,2)))
                                    xt = xticks;
                                    xticklabels(linspace(events(e_i,1),events(e_i,2),length(xt))*parameters.dt)
                                    xlabel('Seconds')
                                    title(strcat('Event #',string(e_i)))
                                end
                                G_str = string(parameters.G_coeff);
                                del_G_str = string(parameters.del_G_sra);
                                I_str = string(parameters.I_strength);
                                title_str = strcat("G_{coeff} =",G_str,"; \Delta_{G_{SRA}} =",del_G_str);
                                sgtitle(title_str)
%                                 savefig(f,strcat(save_path,'sequence_G_',G_str,'_del_',extractBefore(del_G_str,5),'.fig'))
                                try
                                    short_del_G_str = extractBefore(del_G_str,5);
                                catch
                                    short_del_G_str = del_G_str;
                                end
                                saveas(f,strcat(save_path,'/figures/sequence_G_',G_str,'_del_',short_del_G_str,'_I_',I_str,'.jpg'))
                                close(f)
                            %end
                        %end
                    %end
                %end
            end

            %Third value: Average event length
            vec(3) = avg_event_length;
        end
    end    
end