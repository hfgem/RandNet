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
%
% New parallelized parameter grid code, based on simulation function randnet_calculator:
% randnet_param_tests.m 
%   ->parallelize_parameter_tests_2 
%       ->parallelize_networks_2
%           ->parallelize_network_tests_2
%               ->randnet_calculator

%% Save Path + Load Parameters
addpath('functions')

saveFlag = 0 % 1 to save simulation results
selectSavePath = 0; % 1 to select save destination, 0 to save in results dir
selectLoadPath = 0; % 1 to select load source, 0 to load from results dir
plotResults = 1; % 1 to plot basic simulation results
scriptFolder = '/param_tests'; % sub-folder so save analysis results to

% Save path
if selectSavePath
    save_path = uigetdir('/Users/hannahgermaine/Documents/PhD/'); %Have user input where they'd like the output stored
else
    save_path = [pwd, '/results', scriptFolder];
end

if saveFlag & ~isfolder(save_path)
    mkdir(save_path);
end

%Load parameters
% Save path
if selectLoadPath
    load_path = uigetdir('/Users/hannahgermaine/Documents/PhD/');
else
    load_path = [pwd, '/results'];
end

load(strcat(load_path,'/parameters.mat'))
param_names = fieldnames(parameters);


%% Parameters that are different from the loaded parameters
parameters.saveFlag = saveFlag; % needed to override the loaded parameters
parameters.plotResults = plotResults; % needed to override the loaded parameters

parameters.event_cutoff = 0;
parameters.min_avg_fr = 0.001;
parameters.max_avg_fr= inf;
parameters.min_avg_length = 0;
parameters.max_avg_length = inf;

parameters.t_max = 10;
parameters.t_max = 2;

%%
%___________________________________
%____Define dependent parameters____
%___________________________________
parameters.cluster_n = min(parameters.n*2/parameters.clusters,parameters.n); %number of neurons in a cluster (for small n round(n/3), for large n round(n/5)) 

%Interaction constants
parameters.t_steps = parameters.t_max/parameters.dt; %number of timesteps in simulation
parameters.syn_E = parameters.V_syn_E*ones(parameters.n,1); %vector of the synaptic reversal potential for excitatory connections
parameters.syn_I = parameters.V_syn_I*ones(parameters.n,1); %vector of the synaptic reversal potential for inhibitory connections
parameters.IES = ceil(parameters.IEI/parameters.dt); %inter-event-steps = the number of steps to elapse between spikes

%Calculate connection probabilites
parameters.npairs = parameters.n*(parameters.n-1); %total number of possible neuron connections
parameters.nclusterpairs = parameters.cluster_n*(parameters.cluster_n - 1)*parameters.clusters; %total number of possible intra-cluster connections
parameters.cluster_prob = min(parameters.conn_prob*parameters.npairs/parameters.nclusterpairs,1); %0.2041; %intra-cluster connection probability
parameters.n_I = round((1-parameters.p_E)*parameters.n); %number of inhibitory neurons

if saveFlag
    save(strcat(save_path,'/parameters.mat'),'parameters')
end

%% Set Up Grid Search Parameters

%Test parameters
num_nets = 3;
num_inits = 1;

%Number of parameters to test (each)
test_n = 10;
num_params = 3;

% % temp, for testing code
num_nets = 2;
num_inits = 1;
test_n = 5;
% % 

assert(parameters.usePoisson==1)

%Parameter 1: coefficient of input conductance
W_gin_Min = 4.4*10^-9; 
W_gin_Max = 6.4*10^-9; 
W_gin_n = test_n;
W_gin_vec = linspace(W_gin_Min, W_gin_Max, W_gin_n);

%Parameter 2: E-E strength
del_G_syn_E_E_Min = 8.5*10^(-9); 
del_G_syn_E_E_Max = 10.5*10^(-9); 
del_G_syn_E_E_n = test_n;
del_G_syn_E_E_vec = linspace(del_G_syn_E_E_Min, del_G_syn_E_E_Max, del_G_syn_E_E_n);

%Parameter 3: I-E strength
del_G_syn_I_E_Min = 1.3300e-08; 
del_G_syn_I_E_Max = 1.3300e-08; 
del_G_syn_I_E_n = 1;
del_G_syn_I_E_vec = linspace(del_G_syn_I_E_Min, del_G_syn_I_E_Max, del_G_syn_I_E_n);

%Combined into one parameter vector to pass
parameter_vec = [W_gin_vec, del_G_syn_E_E_vec, del_G_syn_I_E_vec];
parameterSets_vec = combvec(W_gin_vec, del_G_syn_E_E_vec, del_G_syn_I_E_vec);


%Save parameter values
if saveFlag
    save(strcat(save_path,'/parameter_vec.mat'),'parameter_vec','-v7.3')
end

%Set up storage matrix
%success = zeros(test_n*ones(1,num_params));
success = zeros(W_gin_n, del_G_syn_E_E_n, del_G_syn_I_E_n);


%% Run Grid Search With Spike Stats Returned


D = parallel.pool.DataQueue;
h = waitbar(0, 'Starting simulation ...');
num_files = size(parameterSets_vec, 2);
% Dummy call to nUpdateWaitbar to initialise
nUpdateWaitbar(num_files, h);
% Go back to simply calling nUpdateWaitbar with the data
afterEach(D, @nUpdateWaitbar);


resultsMat = zeros(size(parameterSets_vec));
resultsStruct = cell(1, size(parameterSets_vec, 2));

tic
parfor ithParamSet = 1:size(parameterSets_vec, 2)
    [resultsMat(:,ithParamSet), resultsStruct{ithParamSet}] = parallelize_parameter_tests_2(...
                parameters, num_nets, num_inits, parameterSets_vec, ithParamSet);
    send(D, 1);
end
runTime = toc

num_spikers = reshape(squeeze(resultsMat(1,:)), W_gin_n, del_G_syn_E_E_n, del_G_syn_I_E_n);
avg_fr = reshape(squeeze(resultsMat(2,:)), W_gin_n, del_G_syn_E_E_n, del_G_syn_I_E_n);
avg_event_length = reshape(squeeze(resultsMat(3,:)), W_gin_n, del_G_syn_E_E_n, del_G_syn_I_E_n);


if saveFlag
    save(strcat(save_path,'/results.mat'),'resultsMat', 'resultsStruct', '-v7.3')
    
    % Save everything, with unique filename based on date-time
    save( strcat(save_path,'/results_', datestr(now,'yyyy-mm-ddTHH-MM'), '.mat'), '-v7.3')
end

%% Visualize Value Grid Search Results

%Recall:
%Parameter 1: coefficient of input conductance (G_coeff)
%Parameter 2: global inhibition strength (I_strength)
%Parameter 3: SRA step size (del_G_sra)

W_gin_vec, del_G_syn_E_E_vec, del_G_syn_I_E_vec

if plotResults
    
    %G_coeff vs I_strength
    num_spikers_G_I = squeeze(mean(num_spikers,3));
    avg_fr_G_I = squeeze(mean(avg_fr,3));
    avg_event_length_G_I = squeeze(mean(avg_event_length,3));
    figure;
    subplot(1,3,1)
    imagesc(W_gin_vec, del_G_syn_E_E_vec, num_spikers_G_I)
    c1 = colorbar(); c1.Label.String = 'Number of Neurons';
    title('Number of Spiking Neurons'); xlabel('W_{EE}'); ylabel('G_{in}')
    subplot(1,3,2)
    imagesc(W_gin_vec, del_G_syn_E_E_vec, avg_fr_G_I)
    c2 = colorbar(); c2.Label.String = "Hz";
    title('Average Firing Rate'); xlabel('W_{EE}'); ylabel('G_{in}')
    subplot(1,3,3)
    imagesc(W_gin_vec, del_G_syn_E_E_vec, avg_event_length_G_I)
    c3 = colorbar(); c3.Label.String = "Seconds";
    title('Average Event Length'); xlabel('W_{EE}'); ylabel('G_{in}')
    clear c1 c2 c3
    
    %G_coeff vs del_G_sra
    num_spikers_G_S = squeeze(mean(num_spikers,2));
    avg_fr_G_S = squeeze(mean(avg_fr,2));
    avg_event_length_G_S = squeeze(mean(avg_event_length,2));
    figure;
    subplot(1,3,1)
    imagesc(W_gin_vec, del_G_syn_I_E_vec, num_spikers_G_S)
    c1 = colorbar(); c1.Label.String = 'Number of Neurons';
    title('Number of Spiking Neurons'); xlabel('W_{IE}'); ylabel('G_{in}')
    subplot(1,3,2)
    imagesc(W_gin_vec, del_G_syn_I_E_vec, avg_fr_G_S)
    c2 = colorbar(); c2.Label.String = "Hz";
    title('Average Firing Rate'); xlabel('W_{IE}'); ylabel('G_{in}')
    subplot(1,3,3)
    imagesc(W_gin_vec, del_G_syn_I_E_vec, avg_event_length_G_S)
    c3 = colorbar(); c3.Label.String = "Seconds";
    title('Average Event Length'); xlabel('W_{IE}'); ylabel('G_{in}')
    clear c1 c2 c3

    %I_strength vs del_G_sra
    num_spikers_I_S = squeeze(mean(num_spikers,1));
    avg_fr_I_S = squeeze(mean(avg_fr,1));
    avg_event_length_I_S = squeeze(mean(avg_event_length,1));
    figure;
    subplot(1,3,1)
    imagesc(del_G_syn_E_E_vec, del_G_syn_I_E_vec, num_spikers_I_S)
    c1 = colorbar(); c1.Label.String = 'Number of Neurons';
    title('Number of Spiking Neurons'); xlabel('W_{IE}'); ylabel('W_{EE}')
    subplot(1,3,2)
    imagesc(del_G_syn_E_E_vec, del_G_syn_I_E_vec, avg_fr_I_S)
    c2 = colorbar(); c2.Label.String = "Hz";
    title('Average Firing Rate'); xlabel('W_{IE}'); ylabel('W_{EE}')
    subplot(1,3,3)
    imagesc(del_G_syn_E_E_vec, del_G_syn_I_E_vec, avg_event_length_I_S)
    c3 = colorbar(); c3.Label.String = "Seconds";
    title('Average Event Length'); xlabel('W_{IE}'); ylabel('W_{EE}')
    clear c1 c2 c3
    
end

%% Functions 


function p = nUpdateWaitbar(data, h)
% https://www.mathworks.com/matlabcentral/answers/660793-help-with-parfor-progress-bar-using-data-queue
persistent TOTAL COUNT H
if nargin == 2
    % initialisation mode
    H = h;
    TOTAL = data;
    COUNT = 0;
else
    % afterEach call, increment COUNT
    COUNT = 1 + COUNT;
    p = COUNT / TOTAL;
    waitbar(p, H, ...
        ['Simulation ', num2str(COUNT), ' of ', num2str(TOTAL), ...
        ' complete. Runtime: ', datestr(datenum(0,0,0,0,0,toc),'HH:MM:SS'), ...
        ' Remaining time: ', datestr(datenum(0,0,0,0,0,toc/p-toc),'HH:MM:SS')]);
end
end
