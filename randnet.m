%RandNet Project - random networks with global inhibition and their ability
%to produce sequences

%% Save Path

clear all

saveFlag = 0; % 1 to save simulation results
selectPath = 0; % 1 to select save destination, 0 to save in current dir
plotResults = 1; % 1 to plot basic simulation results

if selectPath
    save_path = uigetdir('/Users/hannahgermaine/Documents/PhD/','Select Save Folder'); %Have user input where they'd like the output stored
else
    save_path = [pwd, '/results'];
end
addpath('functions')

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
tau_sra = 30*10^(-3); %spike rate adaptation time constant (s)
%If want to have STDP, change connectivity_gain to > 0.
connectivity_gain = 0; %0.005; %amount to increase or decrease connectivity by with each spike (more at the range of 0.002-0.005)

% Input conductance, if using non-Poisson input
G_std = -19*10^-9; % STD of the input conductance G_in, if using randn()
G_mean = 0* 10^-12; % mean of the input conductance G_in, if using randn()

% Poisson input parameters
usePoisson = 1; % 1 to use poisson spike inputs, 0 for randn() input
rG = 500; % input spiking rate, if using poisson inputs
W_gin = 5.4*10^-9; % increase in conductance, if using poisson inputs

%Calculate connection probabilites
conn_prob = 0.08; %set a total desired connection probability
p_E = 0.75; %probability of an excitatory neuron

%Global Inhibition
global_inhib = 1; % if 1, I-cells are not clusterd and have connection probability p_I
p_I = 0.5; % probability of an I cell connecting to any other cell


test_val_max = 1; % How many tests of different initializations to run
include_all = 1; % if a neuron is not in any cluster, take cluster membership from a highly connected neuron

%% Parameters for sequence analysis

IEI = 0.02; %inter-event-interval (s) the elapsed time between spikes to count separate events

bin_width = 5*10^(-3); %5 ms bin

%TEST 1: The number of neurons participating in a sequence must pass a threshold:
event_cutoff = 0.10; %0.25; %fraction of neurons that have to be involved to constitute a successful event

%TEST 2: The firing rate must fall within a realistic range
min_avg_fr = 0.01;
max_avg_fr = 3.0;

% TEST 3: The sequence(s) of firing is(are) within reasonable lengths
min_avg_length = 0.01;
max_avg_length = 0.5;



%% Save parameters to a structure and to computer
w = whos;
parameters = struct;
for a = 1:length(w) 
    parameters.(w(a).name) = eval(w(a).name); 
end
clear w a

%% Create Networks and Check Spike Progression
%Runs through a series of different random number generator seeds to change
%the network connectivity and setup, and then automatically outputs
%sequence data to a folder in save_path

%If uploading a parameter file, uncomment the next line
% load(strcat(save_path,'/parameters.mat'))

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
n_I = round((1-parameters.p_E)*parameters.n); %number of inhibitory neurons
%save for easy calculations
parameters.('npairs') = npairs;
parameters.('nclusterpairs') = nclusterpairs;
parameters.('cluster_prob') = cluster_prob;
parameters.('p_I') = p_I;
parameters.('n_I') = n_I;

if saveFlag & ~isfolder(save_path)
    mkdir(save_path);
end
    
if saveFlag
    save(strcat(save_path,'/parameters.mat'),'parameters'); 
end

%____________________________________________
%___Run Simulations for Different Networks___
%____________________________________________

for i = 1:1%10 %how many different network structures to test
    
    %CREATE NETWORK SAVE PATH
    net_save_path = strcat(save_path,'/network_',string(i));
    if saveFlag & ~isfolder(net_save_path)
        mkdir(net_save_path);
    end
    

    network = create_clusters(parameters, 'seed', i, 'include_all', parameters.include_all, 'global_inhib', parameters.global_inhib);
    if saveFlag
        save(strcat(net_save_path,'/network.mat'),'network');
    end

    
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
        if parameters.usePoisson
            G_in = zeros(parameters.n, parameters.t_steps+1);
            for k = 2:(parameters.t_steps+1)
                G_in(:,k) = G_in(:,k-1)*exp(-parameters.dt/parameters.tau_syn_E);
                G_in(:,k) = G_in(:,k) + parameters.W_gin * [rand(parameters.n, 1) < (parameters.dt*parameters.rG)];
            end
        else
            G_in = (parameters.G_std*randn(parameters.n,parameters.t_steps+1))+parameters.G_mean;
            G_in(G_in<0) = 0;
        end
        parameters.('G_in') = G_in;

        %Create Storage Variables
        V_m = zeros(parameters.n,parameters.t_steps+1); %membrane potential for each neuron at each timestep
        V_m(:,1) = parameters.V_reset + randn([parameters.n,1])*(10^(-3))*sqrt(dt); %set all neurons to baseline reset membrane potential with added noise
        
        seed = j;
        
        %Run model
        [V_m, G_sra, G_syn_E_E, G_syn_I_E, G_syn_E_I, G_syn_I_I, conns] = ...
                randnet_calculator(parameters, seed, network, V_m);
        V_m_var(j).V_m = V_m;
        G_var(j).G_in = G_in;
        G_var(j).G_sra = G_sra;
        G_var(j).G_syn_I_E = G_syn_I_E;
        G_var(j).G_syn_E_E = G_syn_E_E;
        G_var(j).G_syn_I_I = G_syn_I_I;
        G_var(j).G_syn_E_I = G_syn_E_I;
        
        [network_spike_sequences, network_cluster_sequences] = detect_events(parameters, network, V_m , j, network_spike_sequences, network_cluster_sequences);
    
        if plotResults

            %Find spike profile
            spikes_V_m = V_m >= parameters.V_th;
            [spikes_x,spikes_t] = find(spikes_V_m);
            max_time = max(spikes_t);
            spiking_neurons = unique(spikes_x, 'stable');

            if isfield(network_spike_sequences(j), 'events')
                events = network_spike_sequences(j).events;

                %Visualize re-ordered spike sequences
                if any(strcmp('spike_order',fieldnames(network_spike_sequences)))
                    f = figure;
                    axes = [];
                    num_events = size(events, 1)
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
                    sgtitle('Spiking Behavior','FontSize',16)

                    %linkaxes(axes)
                    if saveFlag
                        savefig(f,strcat(net_save_path,'/','_',string(j),'firing_sequence.fig'))
                        saveas(f,strcat(net_save_path,'/', '_',string(j),'firing_sequence.jpg'))
                        close(f)
                    end
                    clear e_i spike_order reordered_spikes event_length s_i ax ...
                        axes xt xtlbl
                end      
            end
            
            t = [0:dt:t_max];
            figure; plot(t, V_m(1:2,:)); ylabel('Vm (V)'); xlabel('Time (s)'); 
            % figure; plot(t, V_m); ylabel('Vm (V)'); xlabel('Time (s)'); 
            figure; plot(t, G_in(1:2,:)); ylabel('G in (S)'); xlabel('Time (s)'); 
            % figure; plot(t, G_in(1:2,:)); ylabel('G in (S)'); xlabel('Time (s)'); 

            figure; hold on
            if exist('events')==1 
                for i = 1:size(events, 1)
                    fill([t(events(i,:)), fliplr(t(events(i,:)))], [0, 0, size(spikes_V_m, 1), size(spikes_V_m, 1)], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
                end
            end
            plotSpikeRaster( spikes_V_m, 'TimePerBin', parameters.dt, 'PlotType', 'scatter'); 
            ylabel('Cell'); xlabel('Time (s)'); 

            movmeanWindow = (1/parameters.dt) * 0.05;
            meanPopRate = movmean(mean(spikes_V_m, 1)/parameters.dt, movmeanWindow);
            figure; hold on
            if exist('events')==1 
                for i = 1:size(events, 1)
                    fill([t(events(i,:)), fliplr(t(events(i,:)))], [0, 0, ceil(max(meanPopRate)), ceil(max(meanPopRate))], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
                end
            end
            plot(t, meanPopRate)
            ylabel('Population mean rate (Hz)'); xlabel('Time (s)'); 
            yline(mean(meanPopRate), 'g')
            yline(mean(meanPopRate)+ std(meanPopRate))
            yline(mean(meanPopRate)+ 2*std(meanPopRate))
        end
        
    end
    
    
    %SAVE NETWORK DATA
    if saveFlag
        save(strcat(net_save_path,'/V_m_var.mat'),'V_m_var','-v7.3')
        save(strcat(net_save_path,'/G_var.mat'),'G_var','-v7.3')
        save(strcat(net_save_path,'/I_var.mat'),'I_var','-v7.3')
        save(strcat(net_save_path,'/network_spike_sequences.mat'),'network_spike_sequences','-v7.3')
        save(strcat(net_save_path,'/network_cluster_sequences.mat'),'network_cluster_sequences','-v7.3')
    end
    
end %End of network structure loop

try; disp(events); end

structfun(@sum, network_spike_sequences.nonspiking_neurons)

