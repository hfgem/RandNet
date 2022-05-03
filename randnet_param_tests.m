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
selectLoadPath = 1; % 1 to select load source, 0 to load from results dir
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
%{
{assert(parameters.E_events_only==1)
parameters.event_cutoff = 0;
parameters.min_avg_fr = 0.001;
parameters.max_avg_fr= inf;
parameters.min_avg_length = 0;
parameters.max_avg_length = inf;
%}

% Simulation duration
%parameters.t_max = 10;
parameters.t_max = 20;

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
num_nets = 10;
% num_nets = 4;
num_inits = 1;
test_n = 25; % Number of parameters to test (each)


% % temp, for testing code
%{
num_nets = 2;
num_inits = 1;
test_n = 4;
%}
% %
assert(parameters.usePoisson==1)


% Parameters must be a field in the parameter structure, and cannot be a
% dependent parameter set in set_depedent_parameters

parameters.W_gin = 750*10^(-12);
%variedParam(1).name = 'W_gin'; % 1st parameter to be varied. Must be a field in the parameter structure
%variedParam(1).range = linspace(550*10^(-12), 950*10^(-12), test_n); % set of values to test param1 at

parameters.del_G_syn_E_E = 750*10^(-12);
%variedParam(2).name = 'del_G_syn_E_E'; % 2nd parameter to be varied
%variedParam(2).range = linspace(550*10^(-12), 950*10^(-12), test_n); % set of values to test param2 at


variedParam(1).name = 'mnc'; % 2nd parameter to be varied
%variedParam(1).range = linspace(1, 21, 81); % set of values to test param2 at
variedParam(1).range = linspace(1, 6, 21); % set of values to test param2 at

variedParam(2).name = 'clusters'; % 2nd parameter to be varied
%variedParam(2).range = [2:1:21]; % set of values to test param2 at
variedParam(2).range = [2:2:42]; % set of values to test param2 at


parameters.del_G_syn_E_I = 500*10^(-12);
parameters.del_G_syn_I_E = 500*10^(-12);

% variedParam(3).name = 'del_G_syn_I_E'; % 2nd parameter to be varied
% variedParam(3).range =  linspace(1.3300e-08, 1.3300e-08, 1); % set of values to test param2 at

% Combine into one parameter vector to pass to parfor function
parameterSets_vec = combvec(variedParam(:).range);


% Exclude cases where mnc>clusters
if isequal(variedParam(1).name, 'mnc') && isequal(variedParam(2).name, 'clusters')
    parameterSets_vec = parameterSets_vec(:,~[parameterSets_vec(1,:)>parameterSets_vec(2,:)]);
end


%% Run Grid Search With Spike Stats Returned

gcp; % starts parallel pool if not already running

% Waitbar code
D = parallel.pool.DataQueue;
h = waitbar(0, 'Starting simulation ...');
num_files = size(parameterSets_vec, 2);
nUpdateWaitbar(num_files, h); % Dummy call to nUpdateWaitbar to initialise
afterEach(D, @nUpdateWaitbar);

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

resultsMat = nan([cellfun(@length, {variedParam.range}), 4]);
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

frac_partic = resultsMat(nSlices{:},1);
avg_fr = resultsMat(nSlices{:},2);
avg_event_length = resultsMat(nSlices{:},3);
avg_n_events = resultsMat(nSlices{:},4);


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

    % paramPlot1 v paramPlot2
    

    %{
    frac_partic = squeeze(mean(frac_partic,3))';
    avg_fr = squeeze(mean(avg_fr,3))';
    avg_event_length = squeeze(mean(avg_event_length,3))';
    avg_n_events = squeeze(mean(avg_n_events,3))';
    %}
    % avg_n_events = resultsMat(nSlices{:},4);
    
    figure;
    subplot(2,2,1)
    imagesc(variedParam(paramPlot1).range, variedParam(paramPlot2).range, frac_partic', 'AlphaData', ~isnan(frac_partic'))
    set(gca,'YDir','normal')
    c1 = colorbar(); c1.Label.String = 'Fraction of neurons';
    title('Frac. firing (event)'); xlabel(variedParam(paramPlot1).name,'Interpreter','none'); ylabel(variedParam(paramPlot2).name,'Interpreter','none')
    
    subplot(2,2,2)
    imagesc(variedParam(paramPlot1).range, variedParam(paramPlot2).range, avg_fr', 'AlphaData', ~isnan(avg_fr'))
    set(gca,'YDir','normal')
    c2 = colorbar(); c2.Label.String = "Hz";
    title('Mean spike rate (trial)'); xlabel(variedParam(paramPlot1).name,'Interpreter','none'); ylabel(variedParam(paramPlot2).name,'Interpreter','none')
    
    subplot(2,2,3)
    imagesc(variedParam(paramPlot1).range, variedParam(paramPlot2).range, avg_event_length', 'AlphaData', ~isnan(avg_event_length'))
    set(gca,'YDir','normal')
    c3 = colorbar(); c3.Label.String = "Seconds";
    title('Mean Event Length'); xlabel(variedParam(paramPlot1).name,'Interpreter','none'); ylabel(variedParam(paramPlot2).name,'Interpreter','none')
    clear c1 c2 c3
    
    subplot(2,2,4)
    imagesc(variedParam(paramPlot1).range, variedParam(paramPlot2).range, avg_n_events'/parameters.t_max, 'AlphaData', ~isnan(avg_n_events'))
    set(gca,'YDir','normal')
    c3 = colorbar(); c3.Label.String = "nEvents / s";
    title('Mean event frequency'); xlabel(variedParam(paramPlot1).name,'Interpreter','none'); ylabel(variedParam(paramPlot2).name,'Interpreter','none')
    clear c1 c2 c3

    
end

%% Plot directly from resultsStruct
%{
paramPlot1 = 1;
paramPlot2 = 2;
% X = nanmean(arrayfun(@(x)mean(x.results.stdRate), resultsStruct), 3);
X = nanmean(arrayfun(@(x)mean(x.results.frac_participation{1}), resultsStruct), 3);
figure; imagesc(variedParam(paramPlot1).range, variedParam(paramPlot2).range, X', 'AlphaData', ~isnan(X'))
set(gca,'YDir','normal')
c2 = colorbar(); c2.Label.String = "Mean Frac participation";
title('Additional plotting'); xlabel(variedParam(paramPlot1).name,'Interpreter','none'); ylabel(variedParam(paramPlot2).name,'Interpreter','none')
%}

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
