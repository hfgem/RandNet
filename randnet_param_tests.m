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
%           ->randnet_calculator_memOpt

%% Save Path + Load Parameters
addpath('functions')

saveFlag = 1; % 1 to save simulation results
selectSavePath = 1; % 1 to select save destination, 0 to save in results dir
selectLoadPath = 1; % 1 to select load source, 0 to load from results dir
plotResults = 1; % 1 to plot basic simulation results
scriptFolder = '/param_tests'; % sub-folder so save analysis results to

% Save path
if selectSavePath
    save_path = uigetdir('/Users/hannahgermaine/Documents/PhD/'); %Have user input where they'd like the output stored
else
    save_path = [pwd, '/results', scriptFolder]; %#ok<UNRCH>
end

if saveFlag && ~isfolder(save_path)
    mkdir(save_path);
end

% Load parameters
if selectLoadPath
    load_path = uigetdir('/Users/hannahgermaine/Documents/PhD/');
else
    load_path = [pwd, '/results']; %#ok<UNRCH>
end
load(strcat(load_path,'/parameters.mat'))


% __Necessary to override the loaded parameters for save preferences__ %
parameters.saveFlag = saveFlag;
parameters.plotResults = plotResults; 
parameters.save_path = save_path; %In the case you'd like to save from within a function

%% Parameters that are different from the loaded parameters

% Analysis parameters: insert below those parameters you want to ensure are
% set for the parameter tests you'll be doing.

% Network structure parameters
parameters.n = 500; %number of neurons
parameters.clusters = 10; % Number of clusters in the network
parameters.mnc = 2; % mean number of clusters each neuron is a member of

try %Ensure that E_events_only = 1
    assert(parameters.E_events_only == 1)
catch
    parameters.E_events_only = 1; %if 1, only consider E-cells for detect_events
end

%Event analysis type
parameters.eventType = 'Seq';

if strcmp(parameters.eventType,'PBE')
    % Analysis parameters for PBE detection
    parameters.PBE_min_Hz = 0.5; % minimum population mean rate during PBE
    parameters.PBE_zscore = 1.0; % minimum stds above mean rate to detect PBE
    parameters.PBE_min_dur = 30 * (1/1000); % minimum duration of a PBE
    parameters.PBE_window =  10 * (1/1000) *(1/parameters.dt); % width of gaussian kernel used to calculate mean pop activity rate
    parameters.PBE_max_combine = 10 * (1/1000); % Combine adjacent PBEs separaeted by less than this duration
else    
    parameters.E_events_only = 1; % if 1, only consider E-cells for detect_events
    parameters.IEI = 0.1; %inter-event-interval (s) the elapsed time between spikes to count separate events
    parameters.bin_width = 5*10^(-3); %5 ms bin

    %TEST 1: The number of neurons participating in a sequence must pass a threshold:
    parameters.event_cutoff = 0.10; %0.25; %fraction of neurons that have to be involved to constitute a successful event

    %TEST 2: The firing rate must fall within a realistic range
    parameters.min_avg_fr = 0.01;
    parameters.max_avg_fr = 3.0;

    % TEST 3: The sequence(s) of firing is(are) within reasonable lengths
    parameters.min_avg_length = 0.01;
    parameters.max_avg_length = 0.5;
end


%% __set/update dependent parameters__ %%
parameters = set_depedent_parameters(parameters);

if saveFlag
    save(strcat(save_path,'/parameters.mat'),'parameters')
end

%% Set Up Grid Search Parameters
% NOTE: Parameters selected for testing cannot be a dependent parameter set in 
% set_depedent_parameters. See randnet.m "Initialize parameters" section
% for a set of independent parameters.

%Use optimized calculator code or not
optFlag = 1; %If 1 uses randnet_calculator_memOpt.m for calculations, if 0 uses randnet_calculator.m

%Perform analyses only on excitatory neurons or all
parameters.E_events_only = 1; %If 1 only analyses excitatory neuron behavior, if 0 analyses entire population's behavior

%Test parameters
parameters.nTrials = 5; % How many tests of different initializations to run
parameters.nNets = 5; % How many networks to run
test_n = 10; % Number of values to test for each parameters

% assert(parameters.usePoisson==1)

%To set up a parameter for testing, use the following format:
%{
variedParam(1).name = 'W_gin'; % 1st parameter to be varied. Must be a field in the parameter structure
variedParam(1).range = linspace(450*10^(-12), 1050*10^(-12), test_n); % set of values to test param1
variedParam(2).name = 'del_G_syn_E_E'; % 2nd parameter to be varied
variedParam(2).range = linspace(450*10^(-12), 1050*10^(-12), test_n); % set of values to test param2
%}

variedParam(1).name = 'del_G_syn_E_E'; % 2nd parameter to be varied
variedParam(1).range = linspace(1*10^(-12), 100*10^(-12), test_n); % set of values to test param1
variedParam(2).name = 'del_G_syn_I_E'; % 2nd parameter to be varied
variedParam(2).range =  linspace(1*10^(-12), 100*10^(-12), test_n); % set of values to test param2
variedParam(3).name = 'del_G_syn_E_I'; % 3rd parameter to be varied
variedParam(3).range =  linspace(1*10^(-12), 100*10^(-12), test_n); % set of values to test param3
variedParam(4).name = 'p_I'; % 3rd parameter to be varied
variedParam(4).range =  linspace(0, 1, test_n); % set of values to test param3
variedParam(5).name = 'G_std'; %4th parameter to be varied
variedParam(5).range = linspace(1*10^(-12), 100*10^(-12), test_n); % set of values to test param4
%parameters.del_G_syn_I_I = 0;

% Combine into one parameter vector to pass to parfor function
parameterSets_vec = combvec(variedParam(:).range);

% Exclude cases where mnc>clusters
if isequal(variedParam(1).name, 'mnc') && isequal(variedParam(2).name, 'clusters')
    parameterSets_vec = parameterSets_vec(:,~[parameterSets_vec(1,:)>parameterSets_vec(2,:)]);
end

% Save cases and updated parameters
if saveFlag
    save(strcat(save_path,'/parameters.mat'),'parameters')
    save(strcat(save_path,'/variedParam.mat'),'variedParam')
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
parfor ithParamSet = 1:size(parameterSets_vec, 2) %Run through all parameter combinations
    %For each combination run parallelize_parameter_tests_2
    [resultsMatLinear(:,ithParamSet), resultsStructLinear{ithParamSet}] = parallelize_parameter_tests_2(...
                parameters, parameterSets_vec, ithParamSet, variedParam, optFlag);
    send(D, 1);
end
runTime = toc;
fprintf('Program Runtime (s) = %.2f\n',runTime)

%% Format results matrix

resultsMat = nan([cellfun(@length, {variedParam.range}), 4]);
resultsStruct = struct;
for i = 1:size(resultsMatLinear, 2)
    
    structIndices = {};
    for ithParam = 1:size(variedParam, 2)
        structIndices{ithParam} = find(parameterSets_vec(ithParam,i)==variedParam(ithParam).range);
    end
    resultsMat(structIndices{:},:) = resultsMatLinear(:,i);
    
    for j = 1:parameters.nNets
        for k = 1:parameters.nTrials
            resultsStruct(structIndices{:}, j, k).results = resultsStructLinear{i}{j}{k};
        end
    end
end

nSlices = repmat({':'},ndims(resultsMat)-1,1); % get n slices, ":", 

frac_partic = resultsMat(nSlices{:},1);
avg_fr = resultsMat(nSlices{:},2);
avg_event_length = resultsMat(nSlices{:},3);
avg_n_events = resultsMat(nSlices{:},4);


if parameters.saveFlag
    save(strcat(save_path,'/parameterSets_vec.mat'),'parameterSets_vec','-v7.3')
    save(strcat(save_path,'/results.mat'),'resultsMat', 'resultsStruct', '-v7.3')
   
    
    clear D h % Don't save pool.DataQueue or waitbar handle
    clear resultsStructLinear resultsMatLinear % don't save redundant data
    
    % Save everything, with unique filename based on date-time
    save( strcat(save_path,'/results_', datestr(now,'yyyy-mm-ddTHH-MM'), '.mat'), '-v7.3') 
end

%% Visualize Value Grid Search Results

%All pairs of parameters to loop through
pair_param_comb = nchoosek(1:size(variedParam,2),2);

%Create save path if plots are to be saved
if parameters.saveFlag
    fig_save_path = strcat(parameters.save_path,'/grid_pair_results');
    if ~isfolder(fig_save_path)
        mkdir(fig_save_path)
    end
end    

for i = 1:length(pair_param_comb)

    if parameters.plotResults

        % select index of parameters to plot against each other
        % Note: below code is not generalized to arbitrary n of variedParam
        paramPlot1 = pair_param_comb(i,1);
        paramPlot2 = pair_param_comb(i,2);

        notparams = setdiff([1:size(variedParam,2)],[paramPlot1,paramPlot2]);

        % paramPlot1 v paramPlot2 across 4 average results

        f = figure;
        subplot(2,2,1)
        imagesc(variedParam(paramPlot1).range, variedParam(paramPlot2).range, squeeze(mean(frac_partic,notparams))', 'AlphaData', ~isnan(squeeze(mean(frac_partic,notparams))'))
        set(gca,'YDir','normal')
        c1 = colorbar(); c1.Label.String = 'Fraction of neurons';
        title('Frac. firing (event)'); xlabel(variedParam(paramPlot1).name,'Interpreter','none'); ylabel(variedParam(paramPlot2).name,'Interpreter','none')

        subplot(2,2,2)
        imagesc(variedParam(paramPlot1).range, variedParam(paramPlot2).range, squeeze(mean(avg_fr,notparams))', 'AlphaData', ~isnan(squeeze(mean(avg_fr,notparams))'))
        set(gca,'YDir','normal')
        c2 = colorbar(); c2.Label.String = "Hz";
        title('Mean spike rate (trial)'); xlabel(variedParam(paramPlot1).name,'Interpreter','none'); ylabel(variedParam(paramPlot2).name,'Interpreter','none')

        subplot(2,2,3)
        imagesc(variedParam(paramPlot1).range, variedParam(paramPlot2).range, squeeze(mean(avg_event_length,notparams))', 'AlphaData', ~isnan(squeeze(mean(avg_event_length,notparams))'))
        set(gca,'YDir','normal')
        c3 = colorbar(); c3.Label.String = "Seconds";
        title('Mean Event Length'); xlabel(variedParam(paramPlot1).name,'Interpreter','none'); ylabel(variedParam(paramPlot2).name,'Interpreter','none')
        clear c1 c2 c3

        subplot(2,2,4)
        imagesc(variedParam(paramPlot1).range, variedParam(paramPlot2).range, squeeze(mean(avg_n_events,notparams))'/parameters.t_max, 'AlphaData', ~isnan(squeeze(mean(avg_n_events,notparams))'))
        set(gca,'YDir','normal')
        c3 = colorbar(); c3.Label.String = "nEvents / s";
        title('Mean event frequency'); xlabel(variedParam(paramPlot1).name,'Interpreter','none'); ylabel(variedParam(paramPlot2).name,'Interpreter','none')
        clear c1 c2 c3

        f.Units = 'normalized';
        f.Position = [0 0 1 1];

        if parameters.saveFlag
            fig_name = strcat('randnet_param_test_plots_parameters_',variedParam(paramPlot1).name,'_v_',variedParam(paramPlot2).name);
            savefig(f,strcat(fig_save_path,'/',fig_name,'.fig'))
            saveas(f,strcat(fig_save_path,'/',fig_name,'.jpg'))
        end    

    end
    
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

