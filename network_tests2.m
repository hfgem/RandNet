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