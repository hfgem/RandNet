%Testing different network properties for sequence generation

%______ABOUT______:
%This code uses a previously successful set of network parameters by
%loading them (lif_network_postrotation.m must have been run to generate
%parameter files) and then modifies one or more of the parameters for 
%testing.
%__________________

%% Initialization

%Select folder of good parameters to use in tests
load_path = uigetdir('/Users/hannahgermaine/Documents/PhD/','Select network.m Load Folder'); %Have user input where they'd like network structure to be loaded from

%Load data to analyze
load(strcat(load_path,'/parameters.mat'))
param_names = fieldnames(parameters);

%How many tests to run
num_tests = 10;

%____USER INPUT VALUE____
%Select which independent parameter(s) to modify
modify_name = {'I_coeff'};
[~,num_mod] = size(modify_name);

parameters.type = 'current';

%Find which parameter indices are being modified
nonmod = find(strcmp(modify_name{1},param_names) == 0); %nonmodified
if num_mod > 1
    for i = 2:num_mod
        false_vals = strcmp(modify_name{i},param_names);
        nonmod = intersect(nonmod,find(false_vals == 0));
    end
end
mod_ind = setdiff(1:length(param_names)',nonmod);
modify_name = param_names(mod_ind); %update to match order

%Store new mod values
new_mod_vals = zeros(num_mod,num_tests);

for i = 1:length(mod_ind)
    param_name = param_names{mod_ind(i)};
    param_val = parameters.(param_name);
    scale = floor(log10(param_val)); %scale of value - to ensure tests are in proper scale
    new_params = linspace(max(param_val - (num_tests/2)*10^scale,0),param_val + (num_tests/2)*10^scale,num_tests);
    new_mod_vals(i,:) = new_params;
end

clear i param_name param_val scale new_params false_vals

%% Run Tests

save_path = uigetdir('/Users/hannahgermaine/Documents/PhD/','Select Network Save Folder'); %Have user input where they'd like the output stored

best_net = 2; %based on prior data, select which rng initialization for a network is best

for mod_i = 1:num_mod %for each modified variable
    mod_name = modify_name{mod_i};
    for test_i = 1:num_tests %for all test values
        
        parameters_copy = parameters;
        parameters_copy.(mod_name) = new_mod_vals(mod_i,test_i);
        
        display(strcat('Current Value =',string(new_mod_vals(mod_i,test_i))))
        
        %____________________________________
        %___Calculate Dependent Parameters___
        %____________________________________
        cluster_n = min(round(parameters_copy.n*2/parameters_copy.clusters),parameters_copy.n); %number of neurons in a cluster (for small n round(n/3), for large n round(n/5)) 
        parameters_copy.('cluster_n') = cluster_n;
        
        %Interaction constants
        t_steps = parameters_copy.t_max/parameters_copy.dt; %number of timesteps in simulation
        E_syn_E = parameters_copy.V_syn_E*ones(parameters_copy.n,1); %vector of the synaptic reversal potential for excitatory connections
        E_syn_I = parameters_copy.V_syn_I*ones(parameters_copy.n,1); %vector of the synaptic reversal potential for inhibitory connections
        IES = ceil(parameters_copy.IEI/parameters_copy.dt); %inter-event-steps = the number of steps to elapse between spikes
        %save for easy calculations
        parameters_copy.('t_steps') = t_steps;
        parameters_copy.('E_syn_E') = E_syn_E;
        parameters_copy.('E_syn_I') = E_syn_I;
        parameters_copy.('IES') = IES;

        %How many tests of different initializations to run
        if strcmp(parameters_copy.type,'cluster')
            test_val_max = parameters_copy.clusters; %Ensures every cluster is initialized
        else
            test_val_max = 10; %This value can be modified as you see fit
        end
        %save for easy calculations
        parameters_copy.('test_val_max') = test_val_max;

        %Adding an input current to all cells (one of the options must be uncommented)
        x_in = [0:parameters_copy.dt:parameters_copy.t_max];
        % %Rhythmic current input: (uncomment if desired)
        % I_coeff = 0; %5.1*10^(-10); %set to 0 for no input current
        % I_Hz = 1; %frequency of input - 1 Hz or 60 bpm, 4-8 Hz for theta rhythm
        % I_in = I_coeff*(0.5 + 0.5*sin(I_Hz*2*pi*x_in)); %Approximately theta wave input current
        %Noisy input: (uncomment if desired)
        I_in = parameters_copy.I_coeff*randn(parameters_copy.n,parameters_copy.t_steps+1)*parameters_copy.I_scale; %Generally use -0.5-0.5 nA stimulus
        %save for easy calculations
        parameters_copy.('x_in') = x_in;
        parameters_copy.('I_in') = I_in;

        %Calculate connection probabilites
        npairs = parameters_copy.n*(parameters_copy.n-1); %total number of possible neuron connections
        nclusterpairs = parameters_copy.cluster_n*(parameters_copy.cluster_n - 1)*parameters_copy.clusters; %total number of possible intra-cluster connections
        cluster_prob = min(parameters_copy.conn_prob*npairs/nclusterpairs,1); %0.2041; %intra-cluster connection probability
        p_I = 1 - parameters_copy.p_E; %probability of an inhibitory neuron
        n_I = round(p_I*parameters_copy.n); %number of inhibitory neurons
        %save for easy calculations
        parameters_copy.('npairs') = npairs;
        parameters_copy.('nclusterpairs') = nclusterpairs;
        parameters_copy.('cluster_prob') = cluster_prob;
        parameters_copy.('p_I') = p_I;
        parameters_copy.('n_I') = n_I;
        
        rng(best_net) %set random number generator for network structure
        
        %CREATE NETWORK SAVE PATH
        net_save_path = strcat(save_path,'/',mod_name,'/network_',string(best_net),'_mod_',mod_name,'_val_',string(new_mod_vals(mod_i,test_i)));
        if ~isfolder(net_save_path)
            mkdir(net_save_path);
        end
        
        %SET UP NETWORK
        [cluster_mat, conns] = create_clusters(parameters_copy, best_net, 1);
        conns_copy = conns; %just a copy of the connections to maintain for reset runs if there's "plasticity"
        %Randomize excitatory and inhibitory connection strengths based on selected
        %probability.
        all_indices = [1:parameters_copy.n];
        I_indices = datasample(all_indices,parameters_copy.n_I,'Replace',false); %indices of inhibitory neurons
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
        
        parameters_copy.('max_fr') = zeros(1,parameters_copy.test_val_max);
        parameters_copy.('avg_fr') = zeros(1,parameters_copy.test_val_max);
        
        for j = 1:parameters_copy.test_val_max
            seed = j;
            
            %Create Storage Variables
            I_syn = zeros(parameters_copy.n,parameters_copy.t_steps+1); %synaptic current emitted by each neuron at each timestep (A)
            %synaptic conductances for each neuron at each timestep
            G_syn_I = zeros(parameters_copy.n,parameters_copy.t_steps+1); %conductance for presynaptic inhibitory (S)
            G_syn_E = zeros(parameters_copy.n,parameters_copy.t_steps+1); %conductance for presynaptic excitatory (S)
            V_m = zeros(parameters_copy.n,parameters_copy.t_steps+1); %membrane potential for each neuron at each timestep
            V_m(:,1) = parameters_copy.V_reset + randn([parameters_copy.n,1])*(10^(-3)); %set all neurons to baseline reset membrane potential with added noise
            G_sra = zeros(parameters_copy.n,parameters_copy.t_steps+1); %refractory conductance for each neuron at each timestep

            %Run model
            [V_m, G_sra, G_syn_I, G_syn_E, I_syn] = lif_sra_calculator_postrotation(...
            parameters_copy, seed, network, I_syn, G_syn_I, G_syn_E, V_m, G_sra);
            network_var(j).V_m = V_m;
            %network_var(j).G_sra = G_sra;
            %network_var(j).G_syn_I = G_syn_I;
            %network_var(j).G_syn_E = G_syn_E;
            %network_var(j).I_syn = I_syn;
            
            %Find spike profile
            spikes_V_m = V_m >= parameters_copy.V_th;
            [spikes_x,spikes_t] = find(spikes_V_m);
            max_time = max(spikes_t);
            
            %Find maximum firing rate + average maximum firing rates of neurons
            all_fr = sum(spikes_V_m,2)/parameters_copy.t_max;
            max_fr = max(all_fr);
            avg_fr = mean(all_fr);
            
            parameters_copy.max_fr(1,j) = max_fr;
            parameters_copy.avg_fr(1,j) = avg_fr;
            
            %USE FR AS AN INDICATOR OF GOOD SEQUENCES
            %In hippocampus, Jadhav lab uses 3 Hz as a place-cell cutoff
            if max_fr >= 1 %and(max_fr <= 3, max_fr >= 1) && and(avg_fr <= 2, avg_fr >= 1/(parameters_copy.t_max+1))
                %Find event times
                events = []; 
                last_start = spikes_t(1);
                last_time = spikes_t(1);
                spike_count = 1;
                for t_i = 2:length(spikes_t)
                    s_i = spikes_t(t_i);
                    if s_i - last_time <= parameters_copy.IES
                        last_time = s_i;
                        spike_count = spike_count + 1;
                    else
                        if (last_start ~= last_time) && (spike_count > parameters_copy.event_cutoff*parameters_copy.n) %weed out events w/ too few spikes
                            events = [events; [last_start, last_time]]; %#ok<AGROW> %add the last range of spikes to the events vector
                        end
                        last_start = s_i;
                        last_time = s_i;
                        spike_count = 1;
                    end
                end
                if (last_start ~= last_time) && (spike_count > parameters.event_cutoff*parameters.n) %weed out events w/ too few spikes
                    events = [events; [last_start, last_time]]; %#ok<AGROW> %add the last interval
                end
                [num_events,~] = size(events);
                %save to both structures
                network_spike_sequences(j).events = events;
                network_cluster_sequences(j).events = events;

                %Find spike sequences
                for e_i = 1:num_events
                    %store spike orders for each event
                    event_spikes = spikes_V_m(:,events(e_i,1):events(e_i,2));
                    [e_spikes_x, ~] = find(event_spikes);
                    spike_order = unique(e_spikes_x,'stable');
                    network_spike_sequences(j).spike_order.(strcat('sequence_',string(e_i))) = spike_order;
                    %store ranks for each neuron
                    ranks_vec = zeros(1,parameters_copy.n);
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
                    reordered_spikes = zeros(size(sub_spikes_V_m));
                    [~,event_length] = size(reordered_spikes);
                    for s_i = 1:length(spike_order)
                        reordered_spikes(s_i,:) = sub_spikes_V_m(spike_order(s_i),:);
                    end 
                    ax = subplot(1,num_events,e_i);
                    axes = [axes, ax];
                    imagesc(reordered_spikes)
                    xticks(round(linspace(1,event_length,20),2)) %20 ticks will be displayed
                    xt = get(gca,'XTick');
                    xtlbl = round(linspace(events(e_i,1)*parameters_copy.dt,events(e_i,2)*parameters_copy.dt,numel(xt)),2);
                    colormap(flip(gray))
                    xlabel('Time (s)','FontSize',16)
                    ylabel('Reordered Neuron Number','FontSize',16)
                    title(strcat('Event #',string(e_i)))
                    set(gca, 'XTick',xt, 'XTickLabel',xtlbl)
                end
                if strcmp(parameters_copy.type,'cluster')
                    sgtitle(strcat('Spiking Behavior: cluster = ',string(j)),'FontSize',16)
                else
                    sgtitle('Spiking Behavior','FontSize',16)
                end
                %linkaxes(axes)
                savefig(f,strcat(net_save_path,'/',parameters_copy.type,'_',string(j),'firing_sequence.fig'))
                saveas(f,strcat(net_save_path,'/',parameters_copy.type,'_',string(j),'firing_sequence.jpg'))
                close(f)
                clear e_i spike_order reordered_spikes event_length s_i ax ...
                    axes xt xtlbl

                %Find cluster sequence per event by moving bin
                bin_width = 5*10^(-3); %5 ms bin
                bin_size = ceil(bin_width/parameters_copy.dt); %number of timesteps to use in a bin
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
            end
        end
        %SAVE NETWORK DATA
        save(strcat(net_save_path,'/parameters.mat'),'parameters_copy')
        save(strcat(net_save_path,'/network_var.mat'),'network_var','-v7.3')
        save(strcat(net_save_path,'/network_spike_sequences.mat'),'network_spike_sequences','-v7.3')
        save(strcat(net_save_path,'/network_cluster_sequences.mat'),'network_cluster_sequences','-v7.3')
    end
end    
