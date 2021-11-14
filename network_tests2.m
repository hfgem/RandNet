%Testing different network properties for sequence generation

%______ABOUT______:
%This code uses a previously successful set of network parameters by
%loading them (lif_network_postrotation.m must have been run to generate
%parameter files) and then grid searches the parameter space of 4
%parameters:
%   1. tau_sra
%   2. del_G_sra
%   3. del_G_syn_E (through del_G_syn_E_E = del_G_syn_E_I
%   4. del_G_syn_I (through del_G_syn_I_I = del_G_syn_I_E
%
%Each parameter set will be tested for 10 networ initializations and 10
%randomized neuron initializations per network for a total of 100 tests of
%each set. The results of each initialization will be compared against
%criteria of a successful output and a score out of 100 will be calculated.
%This score will be stored in a 4D matrix of the parameter space.
%
%The results of the grid search will be visualized and characterized for a
%parameter set description of best results. Depending on the results,
%another grid search can be performed narrowing ranges of parameters
%manually to determine the best parameter space for solutions.
%__________________

%% Note to self

%Test adding global connectivity to inhibitory neurons - up to ~50% of E
%cells perhaps?

%% Save Path + Load Parameters
%Load parameters
msgbox('Select folder where parameters are stored.')
load_path = uigetdir('/Users/hannahgermaine/Documents/PhD/');
load(strcat(load_path,'/parameters.mat'))
param_names = fieldnames(parameters);

%Create save path
msgbox('Select folder to save results.')
save_path = uigetdir('/Users/hannahgermaine/Documents/PhD/');

%___________________________________
%____Define dependent parameters____
%___________________________________
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

%Adding an input conductance to all cells (one of the options must be uncommented)
x_in = [0:parameters.dt:parameters.t_max];
G_in = parameters.G_coeff*randn(parameters.n,parameters.t_steps+1)*parameters.G_scale;
parameters.('x_in') = x_in;
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

%% Set Up Grid Search Parameters
%Number of parameter values to test
test_n = 5;

%Define parameter vectors
tau_sra_vec = linspace(10*10^(-3),60*10^(-3),test_n);
del_G_sra_vec = linspace(1*10^(-9),200*10^(-9),test_n);
del_G_syn_E_vec = linspace(2*10^(-9),15*10^(-9),test_n);
del_G_syn_I_vec = linspace(2*10^(-9),15*10^(-9),test_n);

parameter_vec = [tau_sra_vec; del_G_sra_vec; del_G_syn_E_vec; del_G_syn_I_vec];
clear tau_sra_vec del_G_sra_vec del_G_syn_E_vec del_G_syn_I_vec

%Save parameter values
save(strcat(save_path,'/parameter_vec.mat'),'parameter_vec','-v7.3')

%Set up storage matrix
success = zeros(test_n,test_n,test_n,test_n);
%dim 1 = tau_sra
%dim 2 = del_G_sra
%dim 3 = del_G_syn_E
%dim 4 = del_G_syn_I

%Test parameters
num_nets = 5; %number of network structure initializations
num_inits = 5; %number of firing initializations

%% Run Grid Search
%Start parallel pool for parallelizing the grid search
% parpool(4)

%Loop through all parameter sets
parfevalOnAll(gcp(), @warning, 0, 'off', 'MATLAB:singularMatrix');
parfor ind = 1:test_n^4, success(ind) = parallelize_parameter_tests(parameters,num_nets,...
    num_inits, parameter_vec, test_n, ind); end
save(strcat(save_path,'/success.mat'),'success','-v7.3')

%% Visualize Grid Search Results
%Recall dimensions:
%dim 1 = tau_sra
%dim 2 = del_G_sra
%dim 3 = del_G_syn_E
%dim 4 = del_G_syn_I

%Calculate and plot average values as function of parameter values to find
%low dimensional trends
avg_tau_sra = sum(success,1)/(test_n^3);
avg_del_G_sra = sum(success,2)/(test_n^3);
avg_del_G_syn_E = sum(success,3)/(test_n^3);
avg_del_G_syn_I = sum(success,4)/(test_n^3);

figure;
ax1 = subplot(2,2,1);
plot(parameter_vec(1,:),avg_tau_sra)
title('Avg Tau SRA')
ax2 = subplot(2,2,2);
plot(parameter_vec(2,:),avg_del_G_sra)
title('Avg \del G SRA')
ax3 = subplot(2,2,3);
plot(parameter_vec(3,:),avg_del_G_syn_E)
title('Avg \del G syn_E')
ax4 = subplot(2,2,4);
plot(parameter_vec(4,:),avg_del_G_syn_I)
linkaxes([ax1,ax2,ax3,ax4])

%Calculate and plot 2D results for parameter pairings


%%

% %Loop through tau_sra values
% for tau_sra = 1:test_n
%     parameters.tau_sra = tau_sra_vec(tau_sra);
%     %display(strcat('Current tau_sra index = ',string(tau_sra)))
%     %Loop through del_G_sra values
%     for del_G_sra = 1:test_n
%         parameters.del_G_sra = del_G_sra_vec(del_G_sra);
%         %display(strcat('Current del_G_sra index = ',string(del_G_sra)))
%         %Loop through del_G_syn_E values
%         for del_G_syn_E = 1:test_n
%             parameters.del_G_syn_E_E = del_G_syn_E_vec(del_G_syn_E);
%             parameters.del_G_syn_E_I = del_G_syn_E_vec(del_G_syn_E);
%             %display(strcat('Current del_G_syn_E index = ',string(del_G_syn_E)))
%             %Loop through del_G_syn_I values
%             for del_G_syn_I = 1:test_n
%                 parameters.del_G_syn_I_I = del_G_syn_I_vec(del_G_syn_I);
%                 parameters.del_G_syn_I_E = del_G_syn_I_vec(del_G_syn_I);
%                 disp([tau_sra, del_G_sra, del_G_syn_E, del_G_syn_I])
%                 %tracker for number of successful trials
% %                 passing_trials = 0;
%                 passing_trials = zeros(1,num_nets);
%                 parfevalOnAll(gcp(), @warning, 0, 'off', 'MATLAB:singularMatrix');
%                 parfor i = 1:num_nets, passing_trials(i) = parallelize_networks(parameters,i, num_inits); end
% %                 for i = 1:num_nets
% %                     rng(i) %set random number generator for network structure
% %                     [cluster_mat, conns] = create_clusters(parameters, i, 1);
% %                     all_indices = [1:parameters.n];
% %                     I_indices = datasample(all_indices,parameters.n_I,'Replace',false); %indices of inhibitory neurons
% %                     E_indices = find(~ismember(all_indices,I_indices)); %indices of excitatory neurons
% %                     clear all_indices
% %                     network = struct;
% %                     network(1).cluster_mat = cluster_mat;
% %                     network(1).conns = conns;
% %                     network(1).I_indices = I_indices;
% %                     network(1).E_indices = E_indices;
% %                     clear cluster_mat conns I_indices E_indices
% %                     pass_vec = zeros(1,num_inits);
% %                     %parfor here in 1 min of running goes through 5 parameter sets
% %                     parfevalOnAll(gcp(), @warning, 0, 'off', 'MATLAB:singularMatrix');
% %                     parfor j = 1:num_inits, pass_vec(j) = parallelize_network_tests(parameters,network,j); end
% %                     passing_trials = passing_trials + sum(pass_vec,'all');
% %                     delete(poolobj)
% %                     for j = 1:num_inits
% %                         seed = j;
% %                         %Create Storage Variables
% %                         I_syn = zeros(parameters.n,parameters.t_steps+1); %synaptic current emitted by each neuron at each timestep (A)
% %                         %synaptic conductances for each neuron at each timestep
% %                         G_syn_I = zeros(parameters.n,parameters.t_steps+1); %conductance for presynaptic inhibitory (S)
% %                         G_syn_E = zeros(parameters.n,parameters.t_steps+1); %conductance for presynaptic excitatory (S)
% %                         V_m = zeros(parameters.n,parameters.t_steps+1); %membrane potential for each neuron at each timestep
% %                         V_m(:,1) = parameters.V_reset + randn([parameters.n,1])*parameters.V_m_noise; %set all neurons to baseline reset membrane potential with added noise
% %                         G_sra = zeros(parameters.n,parameters.t_steps+1); %refractory conductance for each neuron at each timestep
% %                         %Run model
% %                         [V_m, ~, ~, ~, ~] = lif_sra_calculator_postrotation(...
% %                             parameters, seed, network, I_syn, G_syn_I, G_syn_E, V_m, G_sra);
% %                         clear I_syn G_syn_I G_syn_E
% %                         
% %                         %Find spike profile
% %                         spikes_V_m = V_m >= parameters.V_th;
% %                         [spikes_x,spikes_t] = find(spikes_V_m);
% %                         spiking_neurons = unique(spikes_x, 'stable');
% %                         passed_barriers = zeros(1,4);
% %                         %First check if the number of spiking neurons fits criteria
% %                         if length(spiking_neurons) >= parameters.event_cutoff*parameters.n
% %                             passed_barriers(1) = 1;
% %                             %Find maximum firing rate + average maximum firing rates of neurons
% %                             all_fr = sum(spikes_V_m(spiking_neurons,:),2)/parameters.t_max; %only look at neurons that fired at least 1x
% %                             max_fr = max(all_fr);
% %                             avg_fr = mean(all_fr);
% % 
% %                             %Test of success (passing_trials = passing_trials + 1;)
% %                             %According to Hirase et al. 2001, neurons fire
% %                             %around 1 Hz on average and have a max around 1.5 Hz
% %                             if and(avg_fr>= 0.02, avg_fr <= 1.5)
% %                                 passed_barriers(2) = 1;
% %                                 %Find firing event times
% %                                 events = [];
% %                                 event_lengths = [];
% %                                 last_start = spikes_t(1);
% %                                 last_time = spikes_t(1);
% %                                 spike_count = 1;
% %                                 for t_i = 2:length(spikes_t)
% %                                     s_i = spikes_t(t_i);
% %                                     if s_i - last_time <= parameters.IES
% %                                         last_time = s_i;
% %                                         spike_count = spike_count + 1;
% %                                     else
% %                                         if (last_start ~= last_time) && (spike_count > parameters.event_cutoff*parameters.n) %weed out events w/ too few spikes
% %                                             events = [events; [last_start, last_time]]; %#ok<AGROW> %add the last range of spikes to the events vector
% %                                             event_lengths(end+1) = (last_time - last_start)*parameters.dt;
% %                                         end
% %                                         last_start = s_i;
% %                                         last_time = s_i;
% %                                         spike_count = 1;
% %                                     end
% %                                 end
% %                                 clear t_i s_i
% %                                 if (last_start ~= last_time) && (spike_count > parameters.event_cutoff*parameters.n) %weed out events w/ too few spikes
% %                                     events = [events; [last_start, last_time]]; %#ok<AGROW> %add the last interval
% %                                     event_lengths(end+1) = (last_time - last_start)*parameters.dt; %#ok<*SAGROW>
% %                                 end
% %                                 avg_event_length = mean(event_lengths);
% %                                 %Ylinen et al 1995 - average SWR event lasts between 50 ms and 120 ms
% %                                 if and(avg_event_length >= 0.05, avg_event_length <= 0.12)
% %                                     passed_barriers(3) = 1;
% %                                     disp('Passed Event Length Cutoff')
% %                                     %Check membrane potentials
% %                                     %check first event for membrane
% %                                     %potential (good enough!)
% %                                     event_V_m = V_m(spiking_neurons,events(1,1):events(1,2));
% %                                     std_V_m = std(event_V_m,1,'all');
% %                                     %Ylinen et al. 1995 1-5 mV V_m oscillation during SWR
% %                                     if and(std_V_m >= 1*10^(-3), std_V_m <= 5*10^(-3))
% %                                         passed_barriers(4) = 1;
% %                                         passing_trials = passing_trials + 1; %Success! Passed all strict criteria!
% %                                     end %membrane potential if loop  
% %                                 end %SWR length if loop
% %                             end %firing rate if loop
% %                         end %number of neurons if loop
% %                         clear V_m spikes_V_m spikes_x spikes_t spiking_neurons all_fr ...
% %                             max_fr avg_fr events event_lengths last_start ...
% %                             last_time spike_count avg_event_length event_V_m std_V_m
% %                     end %loop of initializations
% %                     clear network j
% %                 end %Loop of network structures
%                 passing_trials = sum(passing_trials,'all');
%                 success(tau_sra,del_G_sra,del_G_syn_E,del_G_syn_I) = ...
%                     passing_trials/(num_nets*num_inits); %store fraction of successful trials
%                 display(strcat('Success = ',string(passing_trials/(num_nets*num_inits))))
%                 %Save current progress
%                 save(strcat(save_path,'/success.mat'),'success','-v7.3')
%             end %del_G_syn_I loop
%         end %del_G_syn_E loop
%     end %del_G_sra loop
% end %tau_sra loop

