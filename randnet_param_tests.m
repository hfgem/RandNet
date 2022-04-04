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
%           ->parallelize_network_tests_2
%               ->randnet_calculator

%% Save Path + Load Parameters
addpath('functions')

saveFlag = 1 % 1 to save simulation results
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

% Load parameters
if selectLoadPath
    load_path = uigetdir('/Users/hannahgermaine/Documents/PhD/');
else
    load_path = [pwd, '/results'];
end
load(strcat(load_path,'/parameters.mat'))


%% Parameters that are different from the loaded parameters

% Analysis parameters
assert(parameters.E_events_only==1)
parameters.event_cutoff = 0;
parameters.min_avg_fr = 0.001;
parameters.max_avg_fr= inf;
parameters.min_avg_length = 0;
parameters.max_avg_length = inf;

% Simulation duration
parameters.t_max = 10;
%parameters.t_max = 2;

% __Necessary to override the loaded parameters__ %
parameters.saveFlag = saveFlag;
parameters.plotResults = plotResults; 


%% __set/update dependent parameters__ %%
parameters = set_depedent_parameters(parameters);

if saveFlag
    save(strcat(save_path,'/parameters.mat'),'parameters')
end


%% Set Up Grid Search Parameters

%Test parameters
num_nets = 5;
num_inits = 1;
test_n = 75; % Number of parameters to test (each)


% % temp, for testing code
% num_nets = 2;
% num_inits = 1;
% test_n = 10;
assert(parameters.usePoisson==1)
% %



% Parameters must be a field in the parameter structure, and cannot be a
% dependent parameter set in set_depedent_parameters

parameters.W_gin = 6.5e-9;
% variedParam(1).name = 'W_gin'; % 1st parameter to be varied. Must be a field in the parameter structure
% variedParam(1).range = linspace(3.4*10^-9, 7.4*10^-9, test_n); % set of values to test param1 at

% parameters.del_G_syn_E_E = 9e-9;
variedParam(2).name = 'del_G_syn_E_E'; % 2nd parameter to be varied
variedParam(2).range = linspace(7.0*10^(-9), 12.0*10^(-9), test_n); % set of values to test param2 at

variedParam(1).name = 'clusters'; % 2nd parameter to be varied
variedParam(1).range =  [2:1:21]; % set of values to test param2 at


parameters.del_G_syn_I_E = 1.3300e-08;
% variedParam(3).name = 'del_G_syn_I_E'; % 2nd parameter to be varied
% variedParam(3).range =  linspace(1.3300e-08, 1.3300e-08, 1); % set of values to test param2 at

% Combine into one parameter vector to pass to parfor function
parameterSets_vec = combvec(variedParam(:).range);




%% Run Grid Search With Spike Stats Returned

% Waitbar code
D = parallel.pool.DataQueue;
h = waitbar(0, 'Starting simulation ...');
num_files = size(parameterSets_vec, 2);
nUpdateWaitbar(num_files, h); % Dummy call to nUpdateWaitbar to initialise
afterEach(D, @nUpdateWaitbar);

gcp; % starts parallel pool if not already running
tic
resultsMatLinear = zeros(4, size(parameterSets_vec, 2));
resultsStructLinear = cell(1, size(parameterSets_vec, 2));
parfor ithParamSet = 1:size(parameterSets_vec, 2)
    
    [resultsMatLinear(:,ithParamSet), resultsStructLinear{ithParamSet}] = parallelize_parameter_tests_2(...
                parameters, num_nets, num_inits, parameterSets_vec, ithParamSet, variedParam);
    send(D, 1);
end
runTime = toc

%% Format results matrix

resultsMat = zeros([cellfun(@length, {variedParam.range}), 4]);
resultsStruct = struct;
for i = 1:size(resultsMatLinear, 2)
    
    structIndices = {};
    for ithParam = 1:size(variedParam, 2)
        structIndices{ithParam} = find(parameterSets_vec(ithParam,i)==variedParam(ithParam).range);
    end
    resultsMat(structIndices{:},:) = resultsMatLinear(:,i);
    
    for j = 1:num_nets
        for k = 1:num_inits
            resultsStruct(structIndices{:}, j, k).results = resultsStructLinear{i}{j}{k};
        end
    end
end

nSlices = repmat({':'},ndims(resultsMat)-1,1); % get n slices, ":", 

num_spikers = resultsMat(nSlices{:},1);
avg_fr = resultsMat(nSlices{:},2);
avg_event_length = resultsMat(nSlices{:},3);
nEvents = resultsMat(nSlices{:},4);


if saveFlag
    save(strcat(save_path,'/parameterSets_vec.mat'),'parameterSets_vec','-v7.3')
    save(strcat(save_path,'/results.mat'),'resultsMat', 'resultsStruct', '-v7.3')
   
    
    clear D h % Don't save pool.DataQueue or waitbar handle
    clear resultsStructLinear resultsMatLinear % don't save redundant data
    
    % Save everything, with unique filename based on date-time
    save( strcat(save_path,'/results_', datestr(now,'yyyy-mm-ddTHH-MM'), '.mat'), '-v7.3') 
end

%% Visualize Value Grid Search Results

if plotResults
    
    % select index of parameters to plot against each other
    % Note: below code is not generalized to arbitrary n of variedParam
    paramPlot1 = 1;
    paramPlot2 = 2;

    % 1 v 2
    num_spikers_G_I = squeeze(mean(num_spikers,3));
    avg_fr_G_I = squeeze(mean(avg_fr,3));
    avg_event_length_G_I = squeeze(mean(avg_event_length,3));
    avg_n_events_G_I = squeeze(mean(nEvents,3));
    figure;
    subplot(2,2,1)
    imagesc(variedParam(paramPlot1).range, variedParam(paramPlot2).range, num_spikers_G_I)
    c1 = colorbar(); c1.Label.String = 'Number of Neurons';
    title('Number of Spiking Neurons'); xlabel(variedParam(paramPlot1).name); ylabel(variedParam(paramPlot2).name)
    subplot(2,2,2)
    imagesc(variedParam(paramPlot1).range, variedParam(paramPlot2).range, avg_fr_G_I, 'AlphaData', ~isnan(avg_fr_G_I))
    c2 = colorbar(); c2.Label.String = "Hz";
    title('Average Firing Rate'); xlabel(variedParam(paramPlot1).name); ylabel(variedParam(paramPlot2).name)
    subplot(2,2,3)
    imagesc(variedParam(paramPlot1).range, variedParam(paramPlot2).range, avg_event_length_G_I, 'AlphaData', ~isnan(avg_event_length_G_I))
    c3 = colorbar(); c3.Label.String = "Seconds";
    title('Average Event Length'); xlabel(variedParam(paramPlot1).name); ylabel(variedParam(paramPlot2).name)
    clear c1 c2 c3
    subplot(2,2,4)
    imagesc(variedParam(paramPlot1).range, variedParam(paramPlot2).range, avg_n_events_G_I, 'AlphaData', ~isnan(avg_n_events_G_I))
    c3 = colorbar(); c3.Label.String = "Event count";
    title('Average number of events'); xlabel(variedParam(paramPlot1).name); ylabel(variedParam(paramPlot2).name)
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
        ['Simulation ', num2str(COUNT), ' of ', num2str(TOTAL), ' complete', newline...
         'Runtime: ', datestr(datenum(0,0,0,0,0,toc),'HH:MM:SS'), ...
        ', Remaining time: ', datestr(datenum(0,0,0,0,0,toc/p-toc),'HH:MM:SS')]);
end
end
