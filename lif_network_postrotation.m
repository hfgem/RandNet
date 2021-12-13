%LIF Network Rotation Project Expanded

%% Save Path

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
tau_syn_E = 10*10^(-3); %AMPA/NMDA synaptic decay time constant (s)
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
% del_G_syn_E = 2*10^(-9); %synaptic conductance step following spike (S)
% del_G_syn_I = 8*10^(-9); %synaptic conductance step following spike (S)
del_G_syn_E_E = 9.5*10^(-9); %synaptic conductance step following spike (S)
del_G_syn_I_I = 1.4*del_G_syn_E_E; %synaptic conductance step following spike (S)
del_G_syn_E_I = 9.5*10^(-9); %synaptic conductance step following spike (S)
del_G_syn_I_E = 1.4*del_G_syn_E_E; %synaptic conductance step following spike (S)

%___________________________
del_G_sra = 330e-09; %spike rate adaptation conductance step following spike %ranges from 1-200 *10^(-9) (S)
tau_sra = (-2.4*10^5)*del_G_sra + 0.11; %spike rate adaptation time constant (s)
%If want to have STDP, change connectivity_gain to > 0.
connectivity_gain = 0; %0.005; %amount to increase or decrease connectivity by with each spike (more at the range of 0.002-0.005)
IEI = 0.1; %inter-event-interval (s) the elapsed time between spikes to count separate events

%How spikes are initiated:
%'cluster' sets a cluster to threshold;
%'current' means spikes depend on an input current of I_in; 
%'neuron' sets a random selection of 2% of neurons to threshold
type = 'neuron'; %'current'; %'cluster';

% %Rhythmic current input: (uncomment if desired)
% I_coeff = 0; %5.1*10^(-10); %set to 0 for no input current
% I_Hz = 1; %frequency of input - 1 Hz or 60 bpm, 4-8 Hz for theta rhythm
% I_in = I_coeff*(0.5 + 0.5*sin(I_Hz*2*pi*x_in)); %Approximately theta wave input current
%Noisy input: (uncomment if desired)
% I_coeff = 2.7; %2.7; %set to 0 for no input current
% I_scale = 1*10^(-9); %sets the scale of the current input
% Input conductance
G_coeff = 0; %-40;
G_scale = 1*10^(-9);

%Calculate connection probabilites
conn_prob = 0.08; %set a total desired connection probability
p_E = 0.75; %probability of an excitatory neuron

%Event Statistics
event_cutoff = 0.25; %fraction of neurons that have to be involved to constitute a successful event

%Save parameters to a structure and to computer
w = whos;
parameters = struct;
for a = 1:length(w) 
    parameters.(w(a).name) = eval(w(a).name); 
end
clear w a
save(strcat(save_path,'/parameters.mat'),'parameters')

%% Create networks and test spike progressions
%Runs through a series of different random number generator seeds to change
%the network connectivity and setup, and then automatically outputs
%sequence data to a folder

%If uploading a parameter file, uncomment the next line
load(strcat(save_path,'/parameters.mat'))

%____________________________________
%___Calculate Dependent Parameters___
%____________________________________
cluster_n = min(parameters.n*2/parameters.clusters,parameters.n); %number of neurons in a cluster (for small n round(n/3), for large n round(n/5)) 
parameters.('cluster_n') = cluster_n;

%Interaction constants
t_steps = parameters.t_max/parameters.dt; %number of timesteps in simulation
E_syn_E = parameters.V_syn_E*ones(parameters.n,1); %vector of the synaptic reversal potential for excitatory connections
E_syn_I = parameters.V_syn_I*ones(parameters.n,1); %vector of the synaptic reversal potential for inhibitory connections
IES = ceil(parameters.IEI/parameters.dt); %inter-event-steps = the number of steps to elapse between spikes
%save for easy calculations
parameters.('t_steps') = t_steps;
parameters.('E_syn_E') = E_syn_E;
parameters.('E_syn_I') = E_syn_I;
parameters.('IES') = IES;

%How many tests of different initializations to run
if strcmp(parameters.type,'cluster')
    test_val_max = parameters.clusters; %Ensures every cluster is initialized
else
    test_val_max = 10; %This value can be modified as you see fit
end
%save for easy calculations
parameters.('test_val_max') = test_val_max;

%Adding an input current to all cells (one of the options must be uncommented)
x_in = [0:parameters.dt:parameters.t_max];
% %Rhythmic current input: (uncomment if desired)
% I_coeff = 0; %5.1*10^(-10); %set to 0 for no input current
% I_Hz = 1; %frequency of input - 1 Hz or 60 bpm, 4-8 Hz for theta rhythm
% I_in = I_coeff*(0.5 + 0.5*sin(I_Hz*2*pi*x_in)); %Approximately theta wave input current
%Noisy input: (uncomment if desired)
%I_in = parameters.I_coeff*randn(parameters.n,parameters.t_steps+1)*parameters.I_scale; %Generally use -0.5-0.5 nA stimulus
G_in = parameters.G_coeff*randn(parameters.n,parameters.t_steps+1)*parameters.G_scale;
%save for easy calculations
parameters.('x_in') = x_in;
%parameters.('I_in') = I_in;
parameters.('G_in') = G_in;

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

for i = 1:10 %how many different network structures to test
    rng(i) %set random number generator for network structure
    
    %CREATE NETWORK SAVE PATH
    net_save_path = strcat(save_path,'/network_',string(i));
    if ~isfolder(net_save_path)
        mkdir(net_save_path);
    end
    
    %SET UP NETWORK
    [cluster_mat, conns] = create_clusters(parameters, i, 1);
    conns_copy = conns; %just a copy of the connections to maintain for reset runs if there's "plasticity"
    %Randomize excitatory and inhibitory connection strengths based on selected
    %probability.
    all_indices = [1:parameters.n];
    I_indices = datasample(all_indices,parameters.n_I,'Replace',false); %indices of inhibitory neurons
    E_indices = find(~ismember(all_indices,I_indices)); %indices of excitatory neurons
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
    network_var = struct;
    network_spike_sequences = struct;
    network_cluster_sequences = struct; 
    
    for j = 1:parameters.test_val_max
        seed = j;
        
        %Create Storage Variables
        I_syn = zeros(parameters.n,parameters.t_steps+1); %synaptic current emitted by each neuron at each timestep (A)
        %synaptic conductances for each neuron at each timestep
        G_syn_I = zeros(parameters.n,parameters.t_steps+1); %conductance for presynaptic inhibitory (S)
        G_syn_E = zeros(parameters.n,parameters.t_steps+1); %conductance for presynaptic excitatory (S)
        V_m = zeros(parameters.n,parameters.t_steps+1); %membrane potential for each neuron at each timestep
        V_m(:,1) = parameters.V_reset + randn([parameters.n,1])*(10^(-3)); %set all neurons to baseline reset membrane potential with added noise
        G_sra = zeros(parameters.n,parameters.t_steps+1); %refractory conductance for each neuron at each timestep
        
        %Run model
        [V_m, G_sra, G_syn_I, G_syn_E, I_syn] = lif_sra_calculator_postrotation(...
            parameters, seed, network, I_syn, G_syn_I, G_syn_E, V_m, G_sra);
        network_var(j).V_m = V_m;
        network_var(j).G_sra = G_sra;
        network_var(j).G_syn_I = G_syn_I;
        network_var(j).G_syn_E = G_syn_E;
        network_var(j).I_syn = I_syn;
        
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
                        normalized_cluster_spikes = cluster_spikes ./ sum(cluster_spikes,1);
                        normalized_cluster_spikes(isnan(normalized_cluster_spikes)) = 0;
                        normalized_cluster_mov_sum = cluster_mov_sum ./ sum(cluster_mov_sum,1);
                        normalized_cluster_mov_sum(isnan(normalized_cluster_mov_sum)) = 0;
                        network_cluster_sequences(j).clusters.(strcat('sequence_',string(e_i))) = cluster_spikes;
                        network_cluster_sequences(j).movsum.(strcat('sequence_',string(e_i))) = cluster_mov_sum;
                        network_cluster_sequences(j).normalized_clusters.(strcat('sequence_',string(e_i))) = normalized_cluster_spikes;
                        network_cluster_sequences(j).normalized_cluster_mov_sum.(strcat('sequence_',string(e_i))) = normalized_cluster_mov_sum;
                    end
                    clear bin_width bin_size cluster_spikes cluster_mov_sum e_i
                end %Sequence length loop
            end %Firing rate loop
        end %Number of neurons in sequence loop
    end
    
    %SAVE NETWORK DATA
    save(strcat(net_save_path,'/network_var.mat'),'network_var','-v7.3')
    save(strcat(net_save_path,'/network_spike_sequences.mat'),'network_spike_sequences','-v7.3')
    save(strcat(net_save_path,'/network_cluster_sequences.mat'),'network_cluster_sequences','-v7.3')
end

%% Calculate trajectory similarity
%This code block does not require any of the earlier code blocks to run. It
%simply requires stored data from prior runs of earlier code blocks.

%Select and load data to analyze
%Most networks from current initializations are good.
net_save_path = uigetdir('/Users/hannahgermaine/Documents/PhD/','Select Network Save Folder'); %Have user input what data they want to analyze
slashes = find(net_save_path == '/');
save_path = net_save_path(1:slashes(end));
load(strcat(save_path,'/parameters.mat'))
load(strcat(net_save_path,'/network_spike_sequences.mat'))

n = parameters.n;
inits = parameters.test_val_max;

%Find which simulated data indices have sequences
viable_inits = [];
for i = 1:inits
    try
        isempty(network_spike_sequences(i).spike_order);
        if ~isempty(network_spike_sequences(i).spike_order)
            viable_inits(end+1) = i; %#ok<SAGROW>
        end
    end
end

if length(viable_inits) ~= 1
    %Create a shuffled dataset based on the viable sequences
    shuffle_n = 100;
    %   First generate shuffled data
    [shuffled_spike_sequences] = generate_shuffled_trajectories(n,...
        shuffle_n, network_spike_sequences, viable_inits);
    shuffled_viable_inits = [1:shuffle_n];

    %=====SPEARMAN'S RANK CORRELATIONS=====

    %_____Real Data Ranks_____
    %ranks = Spearman's rank correlation including nonspiking
    %ranks_mod = Spearman's rank correlation excluding nonspiking
    [ranks, ranks_mod] = calculate_trajectory_similarity_spearmans(n, ...
        viable_inits, network_spike_sequences);
    %   Remove NaN Values
    ranks_mod(isnan(ranks_mod)) = 0;
    %   Turn Rank Matrices Into Vectors of Unique Ranks
    ranks_vec = nonzeros(triu(ranks,1)); %all above diagonal
    ranks_mod_vec = nonzeros(triu(ranks_mod,1)); %all above diagonal

    %_____Shuffled Data Ranks_____
    %   Get shuffled ranks
    [shuffled_ranks, shuffled_ranks_mod] = calculate_trajectory_similarity_spearmans(n, ...
        shuffled_viable_inits, shuffled_spike_sequences);
    shuffled_ranks_vec = nonzeros(triu(shuffled_ranks,1)); %all above diagonal
    shuffled_ranks_mod_vec = nonzeros(triu(shuffled_ranks_mod,1)); %all above diagonal

    %_____Calculate Percentiles_____
    ranks_percentile = comp_percentile(shuffled_ranks_vec,mean(ranks_vec));
    ranks_mod_percentile = comp_percentile(shuffled_ranks_mod_vec,mean(ranks_mod_vec));

    %_____Plot Histograms with Percentiles_____
    f = figure;
    ax1 = subplot(1,2,1);
    histogram(shuffled_ranks_vec,'DisplayName','Shuffled Values')
    hold on
    xline(mean(ranks_vec),'-',strcat('Percentile = ',string(ranks_percentile)),'Color','r','DisplayName','Percentile of Average \newline Real Data Value')
    histogram(ranks_vec,'DisplayName','Real Data Values')
    legend()
    title('Spearmans Rank Correlation Rhos with Nonspiking Neurons at End')
    ax2 = subplot(1,2,2);
    histogram(shuffled_ranks_mod_vec,'DisplayName','Shuffled Values')
    hold on
    xline(mean(ranks_mod_vec),'-',strcat('Percentile = ',string(ranks_mod_percentile)),'Color','r','DisplayName','Percentile of Average \newline Real Data Value')
    histogram(ranks_mod_vec,'DisplayName','Real Data Values')
    legend()
    title('Spearmans Rank Correlation Rhos Excluding Nonspiking Neurons')
    f.Position = [187,387,1112,410];
    savefig(f,strcat(net_save_path,'/','spearmans_rank_percentiles.fig'))
    saveas(f,strcat(net_save_path,'/','spearmans_rank_percentiles.jpg'))
    saveas(f,strcat(net_save_path,'/','spearmans_rank_percentiles.svg'))
    close(f)

    %=====MATCHING INDICES=====

    %_____Real Data MIs_____
    %ranks = Spearman's rank correlation including nonspiking
    %ranks_mod = Spearman's rank correlation excluding nonspiking
    [matching_index, matching_index_mod] = calculate_trajectory_similarity_mi(n, ...
        viable_inits, network_spike_sequences);
    %   Remove NaN Values
    matching_index_mod(isnan(matching_index_mod)) = 0;
    %   Turn Rank Matrices Into Vectors of Unique Ranks
    matching_index_vec = nonzeros(triu(matching_index,1)); %all above diagonal
    matching_index_mod_vec = nonzeros(triu(matching_index_mod,1)); %all above diagonal

    %_____Shuffled Data MIs_____
    %   Get shuffled ranks
    [shuffled_matching_index, shuffled_matching_index_mod] = calculate_trajectory_similarity_spearmans(n, ...
        shuffled_viable_inits, shuffled_spike_sequences);
    shuffled_matching_index_vec = nonzeros(triu(shuffled_matching_index,1)); %all above diagonal
    shuffled_matching_index_mod_vec = nonzeros(triu(shuffled_matching_index_mod,1)); %all above diagonal

    %_____Calculate Percentiles_____
    matching_index_percentile = comp_percentile(shuffled_matching_index_vec,mean(matching_index_vec));
    matching_index_mod_percentile = comp_percentile(shuffled_matching_index_mod_vec,mean(matching_index_mod_vec));

    %_____Plot Histograms with Percentiles_____
    f = figure;
    ax1 = subplot(1,2,1);
    histogram(shuffled_matching_index_vec,'DisplayName','Shuffled Values')
    hold on
    xline(mean(matching_index_vec),'-',strcat('Percentile = ',string(matching_index_percentile)),'Color','r','DisplayName','Percentile of Average \newline Real Data Value')
    histogram(matching_index_vec,'DisplayName','Real Data Values')
    legend()
    title('Matching Indices with Nonspiking Neurons at End')
    ax2 = subplot(1,2,2);
    histogram(shuffled_matching_index_mod_vec,'DisplayName','Shuffled Values')
    hold on
    try %#ok<TRYNC>
        xline(mean(matching_index_mod_vec),'-',strcat('Percentile = ',string(matching_index_mod_percentile)),'Color','r','DisplayName','Percentile of Average \newline Real Data Value')
        histogram(matching_index_mod_vec,'DisplayName','Real Data Values')
    end
    legend()
    title('Matching Indices Excluding Nonspiking Neurons')
    f.Position = [187,387,1112,410];
    savefig(f,strcat(net_save_path,'/','MI_percentiles.fig'))
    saveas(f,strcat(net_save_path,'/','MI_percentiles.jpg'))
    saveas(f,strcat(net_save_path,'/','MI_percentiles.svg'))
    close(f)

else
    disp('Only 1 Sequence')
end

%% Visualize network structure
%Select data to visualize
net_save_path = uigetdir('/Users/hannahgermaine/Documents/PhD/','Select Network Save Folder'); %Have user input where they'd like the output stored
load(strcat(net_save_path,'/network.mat'))
slashes = find(net_save_path == '/');
save_path = net_save_path(1:slashes(end));
load(strcat(save_path,'/parameters.mat'))

cluster_mat = network.cluster_mat;
clusters = parameters.clusters;
n = parameters.n;
conns = network.conns;

%First create indices for each neuron in each cluster - staggering
%locations
for i = 1:clusters
    for j = 1:n
        if i <= clusters/2
            x = i;
            y_base = (clusters/2) - abs((clusters/4) - i);
        else
            x = clusters - i;
            y_base = abs((3*clusters/4) - i);
        end
        cluster_plot_indices(i,j,1) = x + 0.1*randn();
        cluster_plot_indices(i,j,2) = y_base + 0.1*randn();
    end
end

%Next find overlapping neuron scatter positions to draw connecting lines
cluster_plot_connections = [];
for i = 1:n
    clusters_in = find(cluster_mat(:,i));
    if length(clusters_in) > 2
        possible_pairs = nchoosek(clusters_in,2);
        for j = 1:length(possible_pairs)
            index_pair = zeros(2,2);
            index_pair(1,:) = cluster_plot_indices(possible_pairs(j,1),i,:);
            index_pair(2,:) = cluster_plot_indices(possible_pairs(j,2),i,:);
            cluster_plot_connections = cat(3,cluster_plot_connections,index_pair); %#ok<AGROW>
        end    
    elseif length(clusters_in) == 2
        index_pair = zeros(2,2);
        index_pair(1,:) = cluster_plot_indices(clusters_in(1),i,:);
        index_pair(2,:) = cluster_plot_indices(clusters_in(2),i,:);
    end
end
[~, ~, cluster_connections_l] = size(cluster_plot_connections);

%Set colors for each cluster
color_options = jet(clusters);

%Plot
f = figure;
leg_items = [];
leg_names = [];
hold on
for i = 1:clusters
    scat = scatter(cluster_plot_indices(i,:,1),cluster_plot_indices(i,:,2),[],color_options(i,:),'filled');
    leg_items(end+1) = scat;
    leg_names = [leg_names, strcat('Cluster #',string(i))]; %#ok<AGROW>
end
for i = 1:cluster_connections_l
    line(cluster_plot_connections(:,1,i),cluster_plot_connections(:,2,i))
end
legend(leg_items,leg_names)
title('Network Diagram','FontSize',16)
set(gca,'xtick',[],'ytick',[])
savefig(f,strcat(net_save_path,'/','network.fig'))
saveas(f,strcat(net_save_path,'/','network.jpg'))
saveas(f,strcat(net_save_path,'/','network.svg'))

%Plot individual cluster connections
f2 = figure;
hold on
for i = 1:clusters
    within_cluster_connections = [];
    %First find and store intra-cluster connections
    neur_ind = find(cluster_mat(i,:));
    for j = neur_ind %row
        for k = neur_ind %column
            if conns(j,k) == 1
                index_pair = zeros(2,2);
                index_pair(1,:) = cluster_plot_indices(i,j,:);
                index_pair(2,:) = cluster_plot_indices(i,k,:);
                within_cluster_connections = cat(3,within_cluster_connections,index_pair); %#ok<AGROW>
            end    
        end
    end
    [~, ~, cluster_connections_l] = size(within_cluster_connections);
    %Plot
    scat = scatter(cluster_plot_indices(i,:,1),cluster_plot_indices(i,:,2),[],color_options(i,:),'filled');
    for i = 1:cluster_connections_l
        line(within_cluster_connections(:,1,i),within_cluster_connections(:,2,i))
    end
end 
title('Intra-Cluster Connectivity','FontSize',16)
set(gca,'xtick',[],'ytick',[])
savefig(f2,strcat(net_save_path,'/','network_clusters.fig'))

%% Creat Plots of Network Properties
%Select data to visualize
net_save_path = uigetdir('/Users/hannahgermaine/Documents/PhD/','Select Network Save Folder'); %Have user input where they'd like the output stored
load(strcat(net_save_path,'/network.mat'))
slashes = find(net_save_path == '/');
save_path = net_save_path(1:slashes(end));
load(strcat(save_path,'/parameters.mat'))

cluster_mat = network.cluster_mat;
clusters = parameters.clusters;
n = parameters.n;
conns = network.conns;

%Plot Cluster Participation as Imagesc
f = figure;
imagesc(cluster_mat)
colorbar('Ticks',[0,1],'TickLabels',{'No','Yes'})
colormap(summer)
xlabel('Neuron Index')
ylabel('Cluster Index')
title('Cluster Participation')
f.Position = [440,720,866,85];
savefig(f,strcat(net_save_path,'/','cluster_participation_visual.fig'))
saveas(f,strcat(net_save_path,'/','cluster_participation_visual.jpg'))
saveas(f,strcat(net_save_path,'/','cluster_participation_visual.svg'))
%close(f)

%Plot Histogram of number of clusters each neuron participates in
num_clusters = sum(cluster_mat,1);
mean_clusters = mean(num_clusters);
f2 = figure;
histogram(num_clusters,'FaceColor',[1.0000,1.0000,0.4000])
xline(mean_clusters,'label',strcat('Mean = ',string(mean_clusters)),...
    'Color',[0,0.5000,0.4000],'LineWidth',1)
ylabel('Number of Neurons')
xlabel('Number of Clusters')
title('Most Neurons Participate in Multiple Clusters')
savefig(f2,strcat(net_save_path,'/','cluster_participation_histogram.fig'))
saveas(f2,strcat(net_save_path,'/','cluster_participation_histogram.jpg'))
saveas(f2,strcat(net_save_path,'/','cluster_participation_histogram.svg'))
%close(f2)

%Plot Histogram of number of neurons that participate in a cluster overlap
pairs = nchoosek(1:clusters,2);
overlap = zeros(1,length(pairs));
for i = 1:length(pairs)
    overlap(1,i) = sum(sum(cluster_mat(pairs(i,:),:),1) == 2); %neurons in both clusters
end
mean_overlap = mean(overlap);
f3 = figure;
histogram(overlap,'FaceColor',[1.0000,0.7812,0.4975])
xline(mean_overlap,'label',strcat('Mean = ',string(mean_overlap)),...
    'Color',[0,0.5000,0.4000],'LineWidth',1)
ylabel('Number of Cluster Pairs')
xlabel('Number of Neurons')
title('Number of Neurons in Cluster Overlap')
savefig(f3,strcat(net_save_path,'/','cluster_overlap_histogram.fig'))
saveas(f3,strcat(net_save_path,'/','cluster_overlap_histogram.jpg'))
saveas(f3,strcat(net_save_path,'/','cluster_overlap_histogram.svg'))
%close(f3)

clear f f2
%% Creat plot of overlap in 2 clusters

%Select data to visualize
net_save_path = uigetdir('/Users/hannahgermaine/Documents/PhD/','Select Network Save Folder'); %Have user input where they'd like the output stored
load(strcat(net_save_path,'/network.mat'))
slashes = find(net_save_path == '/');
save_path = net_save_path(1:slashes(end));
load(strcat(save_path,'/parameters.mat'))

cluster_mat = network.cluster_mat;
clusters = parameters.clusters;
cluster_n = parameters.cluster_n;
n = parameters.n;
conns = network.conns;
clust_ind = randperm(clusters,2); %Randomly pick 2 clusters for visualization
neur_ind = find(sum(cluster_mat(clust_ind,:),1) == 2); %neurons in both clusters
colors = [[0,0.5000,0.7500];[1.0000,0.7812,0.4975];[0,0.5000,0.4000]];

%Store indices and colors of scatter points
plot_ind = zeros(n,2); %index of each neuron scatter plot
scatter_col = zeros(n,3); %color of each neuron scatter plot

%Assign neuron color and location based on its cluster participation
for i = 1:n
    if sum(neur_ind == i) %is in both clusters
        x_center = 1.5;
        y_center = 1;
        scatter_col(i,:) = colors(3,:);
        plot_ind(i,1) = x_center + 0.2*(rand() - 0.5);
        plot_ind(i,2) = y_center + 0.75*(rand() - 0.5);
    else
        in_clust = find(cluster_mat(clust_ind,i));
        if ~isempty(in_clust) %neuron is in one of the clusters
            if in_clust == 2 %set centers of point clouds
                x_center = 2;
                y_center = 1;
                scatter_col(i,:) = colors(2,:);
            else
                x_center = 1;
                y_center = 1;
                scatter_col(i,:) = colors(1,:);
            end
            plot_ind(i,1) = x_center + 0.5*(rand() - 0.5);
            plot_ind(i,2) = y_center + 0.5*(rand() - 0.5);
        end
    end
end

%Plot Cluster Connectivities
f = figure;
hold on
for i = 1:2
    c_i = clust_ind(i);
    line_color = colors(i,:);
    for j = 1:n
        if cluster_mat(c_i,j) == 1
            for k = 1:n
                if cluster_mat(c_i,k) == 1
                    if logical(conns(j,k)) == 1
                        line([plot_ind(j,1),plot_ind(k,1)],[plot_ind(j,2),...
                            plot_ind(k,2)],'LineStyle',':','Color',line_color)
                    end
                end
            end
        end
    end
end
scatter(plot_ind(:,1),plot_ind(:,2),50,scatter_col,'filled')
set(gca,'xtick',[],'ytick',[])
savefig(f,strcat(net_save_path,'/','cluster_overlap_vis.fig'))
saveas(f,strcat(net_save_path,'/','cluster_overlap_vis.jpg'))
saveas(f,strcat(net_save_path,'/','cluster_overlap_vis.svg'))


