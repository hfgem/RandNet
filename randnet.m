%RandNet Project - random networks with global inhibition and their ability
%to produce sequences

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
t_max = 2; %maximum amount of time (s)
dt = 0.1*10^(-3); %timestep (s)
tau_syn_E = 2*10^(-3); %AMPA synaptic decay time constant (s) [Ignoring NMDA as slow and weak]
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
del_G_syn_E_E = 9.5*10^(-9); %synaptic conductance step following spike (S)
del_G_syn_I_I = 1.4*del_G_syn_E_E; %synaptic conductance step following spike (S)
del_G_syn_E_I = del_G_syn_E_E; %synaptic conductance step following spike (S)
del_G_syn_I_E = 1.4*del_G_syn_E_E; %synaptic conductance step following spike (S)

%___________________________
del_G_sra = 330e-09; %spike rate adaptation conductance step following spike %ranges from 1-200 *10^(-9) (S)
tau_sra = (-2.4*10^5)*del_G_sra + 0.11; %spike rate adaptation time constant (s)
%If want to have STDP, change connectivity_gain to > 0.
connectivity_gain = 0; %0.005; %amount to increase or decrease connectivity by with each spike (more at the range of 0.002-0.005)
IEI = 0.05; %inter-event-interval (s) the elapsed time between spikes to count separate events

% Input conductance
G_coeff = -38; %-40;
G_scale = 1*10^(-9);

%Calculate connection probabilites
conn_prob = 0.08; %set a total desired connection probability
p_E = 0.75; %probability of an excitatory neuron

%Global Inhibition
%This value controls the amount of global inhibition - set to 0 if you want
%no global inhibition, and scale upward for scaled connectivity.
I_strength = 1;
only_global = 0; %Set this flag to 1 if you want global inhibition to override random inhibitory connections given by "create_cluters"

%How many tests of different initializations to run
test_val_max = 1; %This value can be modified as you see fit

%Event Statistics
event_cutoff = 0.25; %0.25; %fraction of neurons that have to be involved to constitute a successful event

%Save parameters to a structure and to computer
w = whos;
parameters = struct;
for a = 1:length(w) 
    parameters.(w(a).name) = eval(w(a).name); 
end
clear w a
save(strcat(save_path,'/parameters.mat'),'parameters')

%% Create Networks and Check Spike Progression
%Runs through a series of different random number generator seeds to change
%the network connectivity and setup, and then automatically outputs
%sequence data to a folder in save_path

%If uploading a parameter file, uncomment the next line
load(strcat(save_path,'/parameters.mat'))

%____________________________________
%___Calculate Dependent Parameters___
%____________________________________
cluster_n = min(parameters.n*2/parameters.clusters,parameters.n); %number of neurons in a cluster (for small n round(n/3), for large n round(n/5)) 
parameters.('cluster_n') = cluster_n;

%Interaction constants
t_steps = parameters.t_max/parameters.dt; %number of timesteps in simulation
syn_E = parameters.V_syn_E*ones(parameters.n,1); %vector of the synaptic reversal potential for excitatory connections
syn_I = parameters.V_syn_I*ones(parameters.n,1); %vector of the synaptic reversal potential for inhibitory connections
IES = ceil(parameters.IEI/parameters.dt); %inter-event-steps = the number of steps to elapse between spikes
%save for easy calculations
parameters.('t_steps') = t_steps;
parameters.('syn_E') = syn_E;
parameters.('syn_I') = syn_I;
parameters.('IES') = IES;

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
    rng(i,'twister') %set random number generator for network structure
    
    %CREATE NETWORK SAVE PATH
    net_save_path = strcat(save_path,'/network_',string(i));
    if ~isfolder(net_save_path)
        mkdir(net_save_path);
    end
    
    %SET UP NETWORK
    %Decide which neurons are inhib and excit 
    all_indices = [1:parameters.n];
    I_indices = datasample(all_indices,parameters.n_I,'Replace',false); %indices of inhibitory neurons
    E_indices = find(~ismember(all_indices,I_indices)); %indices of excitatory neurons
    [cluster_mat, conns] = create_clusters(parameters, i, 1);
    %Add in global inhibition, added to individual connections already
    %given. If global inhibition overrides any pre-set connections with
    %inhibitory neurons, reset values to global inhibition values.
    if parameters.only_global
        conns(I_indices,:) = parameters.I_strength*(rand([parameters.n*(1-parameters.p_E), parameters.n]) < parameters.p_I);
    else
        conns(I_indices,:) = conns(I_indices,:) + parameters.I_strength*(rand([parameters.n*(1-parameters.p_E), parameters.n]) < parameters.p_I);
    end
    conns_copy = conns; %just a copy of the connections to maintain for reset runs if there's "plasticity"
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
    V_m_var = struct;
    G_var = struct;
    I_var = struct;
    network_spike_sequences = struct;
    network_cluster_sequences = struct; 
    
    for j = 1:parameters.test_val_max        
        
        %Create input conductance variable
        G_in = parameters.G_coeff*randn(parameters.n,parameters.t_steps+1)*parameters.G_scale;
        parameters.('G_in') = G_in;
        
        %Create Storage Variables
        V_m = zeros(parameters.n,parameters.t_steps+1); %membrane potential for each neuron at each timestep
        V_m(:,1) = parameters.V_reset + randn([parameters.n,1])*(10^(-3))*sqrt(dt); %set all neurons to baseline reset membrane potential with added noise
        
        seed = j;
        
        %Run model
        [V_m, G_sra, G_syn_E_E, G_syn_I_E, G_syn_E_I, G_syn_I_I, G_in, conns] = ...
                randnet_calculator(parameters, seed, network, V_m);
        V_m_var(j).V_m = V_m;
        G_var(j).G_in = G_in;
        G_var(j).G_sra = G_sra;
        G_var(j).G_syn_I_E = G_syn_I_E;
        G_var(j).G_syn_E_E = G_syn_E_E;
        G_var(j).G_syn_I_I = G_syn_I_I;
        G_var(j).G_syn_E_I = G_syn_E_I;
        
        %Find spike profile
        spikes_V_m = V_m >= parameters.V_th;
        [spikes_x,spikes_t] = find(spikes_V_m);
        max_time = max(spikes_t);
        spiking_neurons = unique(spikes_x, 'stable');
        
        %TEST 1: The number of neurons participating in a sequence must
        %pass a threshold:
        if length(spiking_neurons) >= parameters.event_cutoff*parameters.n
        
            %Find maximum firing rate + average maximum firing rates of neurons
            all_fr = sum(spikes_V_m,2)/parameters.t_max;
            max_fr = max(all_fr);
            avg_fr = mean(all_fr);
            display(avg_fr)

            %TEST 2: The firing rate must fall within a realistic range
            if and(avg_fr>= 0.02, avg_fr <= 1.5)
                %Find event times
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
                            events = [events; [last_start, last_time]]; %#ok<AGROW> %add the last range of spikes to the events vector
                            event_lengths = [event_lengths, (last_time - last_start)*parameters.dt]; %#ok<*AGROW>
                        end
                        last_start = s_i;
                        last_time = s_i;
                        spike_count = 1;
                    end
                end
                if (last_start ~= last_time) && (spike_count > parameters.event_cutoff*parameters.n) %weed out events w/ too few spikes
                    events = [events; [last_start, last_time]]; %#ok<AGROW> %add the last interval
                    event_lengths = [event_lengths, (last_time - last_start)*parameters.dt]; %#ok<*SAGROW>
                end
                avg_event_length = mean(event_lengths);
                [num_events,~] = size(events);
                %save to both structures
                network_spike_sequences(j).events = events;
                network_spike_sequences(j).event_lengths = event_lengths;
                network_cluster_sequences(j).events = events;
                network_cluster_sequences(j).event_lengths = event_lengths;
                
                %TEST 3: The sequence(s) of firing is(are) within
                %reasonable lengths.
                if and(avg_event_length >= 0.02, avg_event_length <= 0.15)

                    %Find spike sequences
                    for e_i = 1:num_events
                        %store spike orders for each event
                        event_spikes = spikes_V_m(:,events(e_i,1):events(e_i,2));
                        [e_spikes_x, ~] = find(event_spikes);
                        spike_order = unique(e_spikes_x,'stable');
                        network_spike_sequences(j).spike_order.(strcat('sequence_',string(e_i))) = spike_order;
                        %store ranks for each neuron
                        ranks_vec = zeros(1,parameters.n);
                        for k = 1:length(spike_order)
                            n_ind = spike_order(k);
                            ranks_vec(1,n_ind) = k;
                        end
                        network_spike_sequences(j).spike_ranks.(strcat('sequence_',string(e_i))) = ranks_vec;
                        %store nonspiking neurons
                        nonspiking_neurons = isnan(ranks_vec./ranks_vec);
                        network_spike_sequences(j).nonspiking_neurons.(strcat('sequence_',string(e_i))) = nonspiking_neurons;
                    end
                    clear e_i event_spikes e_spikes_x spike_order ranks_vec k n_ind nonspiking_neurons

                    %Visualize re-ordered spike sequences
                    f = figure;
                    axes = [];
                    for e_i = 1:num_events
                        spike_order = network_spike_sequences(j).spike_order.(strcat('sequence_',string(e_i)));
                        sub_spikes_V_m = spikes_V_m(:,events(e_i,1):events(e_i,2));
                        reordered_spikes = sub_spikes_V_m(spike_order,:);
                        [~,event_length] = size(reordered_spikes); 
                        ax = subplot(1,num_events,e_i);
                        axes = [axes, ax];
                        imagesc(reordered_spikes)
                        xticks(round(linspace(1,event_length,20))) %20 ticks will be displayed
                        xt = get(gca,'XTick');
                        xtlbl = round(linspace(events(e_i,1)*parameters.dt,events(e_i,2)*parameters.dt,numel(xt)),2);
                        colormap(flip(gray))
                        xlabel('Time (s)','FontSize',16)
                        ylabel('Reordered Neuron Number','FontSize',16)
                        title(strcat('Event #',string(e_i)))
                        set(gca, 'XTick',xt, 'XTickLabel',xtlbl)
                    end
                    if strcmp(parameters.type,'cluster')
                        sgtitle(strcat('Spiking Behavior: cluster = ',string(j)),'FontSize',16)
                    else
                        sgtitle('Spiking Behavior','FontSize',16)
                    end
                    %linkaxes(axes)
                    savefig(f,strcat(net_save_path,'/',parameters.type,'_',string(j),'firing_sequence.fig'))
                    saveas(f,strcat(net_save_path,'/',parameters.type,'_',string(j),'firing_sequence.jpg'))
                    close(f)
                    clear e_i spike_order reordered_spikes event_length s_i ax ...
                        axes xt xtlbl

                    %Find cluster sequence per event by moving bin
                    bin_width = 5*10^(-3); %5 ms bin
                    bin_size = ceil(bin_width/parameters.dt); %number of timesteps to use in a bin
                    for e_i = 1:num_events
                        cluster_spikes = network.cluster_mat*spikes_V_m(:,events(e_i,1):events(e_i,2));
                        cluster_mov_sum = movsum(cluster_spikes',bin_size)';
                        network_cluster_sequences(j).clusters.(strcat('sequence_',string(e_i))) = cluster_spikes;
                        network_cluster_sequences(j).movsum.(strcat('sequence_',string(e_i))) = cluster_mov_sum;
                    end
                    clear bin_width bin_size cluster_spikes cluster_mov_sum e_i
                end %Sequence length loop
            end %Firing rate loop
        end %Number of neurons in sequence loop
    end
    
    %SAVE NETWORK DATA
    save(strcat(net_save_path,'/V_m_var.mat'),'V_m_var','-v7.3')
    save(strcat(net_save_path,'/G_var.mat'),'G_var','-v7.3')
    save(strcat(net_save_path,'/I_var.mat'),'I_var','-v7.3')
    save(strcat(net_save_path,'/network_spike_sequences.mat'),'network_spike_sequences','-v7.3')
    save(strcat(net_save_path,'/network_cluster_sequences.mat'),'network_cluster_sequences','-v7.3')

    
end %End of network structure loop
