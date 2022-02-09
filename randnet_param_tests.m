%Testing different network parameters for sequence generation

%______ABOUT______:
%This code uses a baseline set of network parameters by
%loading them (lif_network_postrotation.m must have been run to generate
%parameter files) and then grid searches the parameter space of defined
%parameters.
%
%Each parameter set will be tested for num_nets network initializations and
%num_inits randomized neuron initializations per network. The results of 
%each initialization will be compared against criteria of a successful
%output and a fractional score will be calculated. This score will be 
%stored in a 4D matrix of the parameter space.
%
%The results of the grid search will be visualized and characterized for a
%parameter set description of best results. Depending on the results,
%another grid search can be performed narrowing ranges of parameters
%manually to determine the best parameter space for solutions.
%__________________

%% Save Path + Load Parameters
%Load parameters
msgbox('Select folder where parameters are stored.')
load_path = uigetdir('/Users/hannahgermaine/Documents/PhD/');
load(strcat(load_path,'/parameters.mat'))
param_names = fieldnames(parameters);

%Create save path
msgbox('Select folder to save results.')
save_path = uigetdir('/Users/hannahgermaine/Documents/PhD/');

%Create figures folder
mkdir(strcat(save_path,'/figures/'))

%___________________________________
%____Define dependent parameters____
%___________________________________
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

%% Set Up Grid Search Parameters

%Test parameters
num_nets = 1;
num_inits = 1;

%Number of parameters to test (each)
test_n = 5;

%Parameter 1: coefficient of input conductance
G_coeff_vec = linspace(-100,100,test_n);

%Parameter 2: global inhibition strength
I_strength_vec = linspace(0,1,test_n);

%Parameter 3: SRA step size
del_G_sra_vec = linspace(1*10^(-9),200*10^(-9),test_n);

%Combined into one parameter vector to pass
parameter_vec = [G_coeff_vec; I_strength_vec; del_G_sra_vec];
clear G_coeff_vec I_strength_vec

%Save parameter values
save(strcat(save_path,'/parameter_vec.mat'),'parameter_vec','-v7.3')

%Set up storage matrix
[num_params, ~] = size(parameter_vec);
success = zeros(test_n*ones(1,num_params));

%% Run Grid Search With Spike Stats Returned
%Start parallel pool for parallelizing the grid search
% parpool(4)

results = zeros(test_n^3,3);

%Loop through all parameter sets
parfevalOnAll(gcp(), @warning, 0, 'off', 'MATLAB:singularMatrix');
parfor ind = 1:test_n^3, results(ind,:) = parallelize_parameter_tests_2(parameters,num_nets,...
    num_inits, parameter_vec, test_n, ind, save_path); end
save(strcat(save_path,'/results.mat'),'results','-v7.3')

num_spikers = reshape(squeeze(results(:,1)),test_n,test_n,test_n);
avg_fr = reshape(squeeze(results(:,2)),test_n,test_n,test_n);
avg_event_length = reshape(squeeze(results(:,3)),test_n,test_n,test_n);

%% Visualize Value Grid Search Results

%Recall:
%Parameter 1: coefficient of input conductance (G_coeff)
%Parameter 2: global inhibition strength (I_strength)
%Parameter 3: SRA step size (del_G_sra)

%G_coeff vs I_strength
num_spikers_G_I = squeeze(mean(num_spikers,3));
avg_fr_G_I = squeeze(mean(avg_fr,3));
avg_event_length_G_I = squeeze(mean(avg_event_length,3));
figure;
subplot(1,3,1)
imagesc(num_spikers_G_I())
c1 = colorbar();
c1.Label.String = 'Number of Neurons';
title('Number of Spiking Neurons')
xticks(1:test_n)
xticklabels(parameter_vec(1,:))
yticks(1:test_n)
yticklabels(parameter_vec(2,:))
xlabel('G_{coeff}')
ylabel('I_{strength}')
subplot(1,3,2)
imagesc(avg_fr_G_I)
c2 = colorbar();
c2.Label.String = "Hz";
title('Average Firing Rate')
xticks(1:test_n)
xticklabels(parameter_vec(1,:))
yticks(1:test_n)
yticklabels(parameter_vec(2,:))
xlabel('G_{coeff}')
ylabel('I_{strength}')
subplot(1,3,3)
imagesc(avg_event_length_G_I)
c3 = colorbar();
c3.Label.String = "Seconds";
title('Average Event Length')
xticks(1:test_n)
xticklabels(parameter_vec(1,:))
yticks(1:test_n)
yticklabels(parameter_vec(2,:))
xlabel('G_{coeff}')
ylabel('I_{strength}')

clear c1 c2 c3

%G_coeff vs del_G_sra
num_spikers_G_S = squeeze(mean(num_spikers,2));
avg_fr_G_S = squeeze(mean(avg_fr,2));
avg_event_length_G_S = squeeze(mean(avg_event_length,2));
figure;
subplot(1,3,1)
imagesc(num_spikers_G_S)
c1 = colorbar();
c1.Label.String = 'Number of Neurons';
title('Number of Spiking Neurons')
xticks(1:test_n)
xticklabels(parameter_vec(1,:))
yticks(1:test_n)
yticklabels(parameter_vec(3,:))
xlabel('G_{coeff}')
ylabel('\delta G_{SRA}')
subplot(1,3,2)
imagesc(avg_fr_G_S)
c2 = colorbar();
c2.Label.String = "Hz";
title('Average Firing Rate')
xticks(1:test_n)
xticklabels(parameter_vec(1,:))
yticks(1:test_n)
yticklabels(parameter_vec(3,:))
xlabel('G_{coeff}')
ylabel('\delta G_{SRA}')
subplot(1,3,3)
imagesc(avg_event_length_G_S)
c3 = colorbar();
c3.Label.String = "Seconds";
title('Average Event Length')
xticks(1:test_n)
xticklabels(parameter_vec(1,:))
yticks(1:test_n)
yticklabels(parameter_vec(3,:))
xlabel('G_{coeff}')
ylabel('\delta G_{SRA}')

clear c1 c2 c3

%I_strength vs del_G_sra
num_spikers_I_S = squeeze(mean(num_spikers,1));
avg_fr_I_S = squeeze(mean(avg_fr,1));
avg_event_length_I_S = squeeze(mean(avg_event_length,1));
figure;
subplot(1,3,1)
imagesc(num_spikers_I_S)
c1 = colorbar();
c1.Label.String = 'Number of Neurons';
title('Number of Spiking Neurons')
xticks(1:test_n)
xticklabels(parameter_vec(2,:))
yticks(1:test_n)
yticklabels(parameter_vec(3,:))
xlabel('I_{strength}')
ylabel('\delta G_{SRA}')
subplot(1,3,2)
imagesc(avg_fr_I_S)
c2 = colorbar();
c2.Label.String = "Hz";
title('Average Firing Rate')
xticks(1:test_n)
xticklabels(parameter_vec(2,:))
yticks(1:test_n)
yticklabels(parameter_vec(3,:))
xlabel('I_{strength}')
ylabel('\delta G_{SRA}')
subplot(1,3,3)
imagesc(avg_event_length_I_S)
c3 = colorbar();
c3.Label.String = "Seconds";
title('Average Event Length')
xticks(1:test_n)
xticklabels(parameter_vec(2,:))
yticks(1:test_n)
yticklabels(parameter_vec(3,:))
xlabel('I_{strength}')
ylabel('\delta G_{SRA}')

clear c1 c2 c3