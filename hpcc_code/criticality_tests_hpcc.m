%HPCC criticality test code
%Pulled from randnet_param_tests, this code is written for the Brandeis
%High Performance Computing Cluster (HPCC). The focus of this code is to
%test different parameter combinations for network criticality.

my_dir = '/home/hgerm/code';
cd(my_dir)
addpath(genpath(my_dir))

%Set Save Parameters
saveFlag = 1; % 1 to save simulation results
selectSavePath = 1; % 1 to select save destination, 0 to save in results dir
selectLoadPath = 1; % 1 to select load source, 0 to load from results dir
plotResults = 1; % 1 to plot basic simulation results
save_path = [my_dir, '/results'];
if saveFlag && ~isfolder(save_path) %Create save directory if doesn't exist
    mkdir(save_path);
end

% Load Parameters
load(strcat(my_dir,'/parameters.mat')) 
parameters.saveFlag = saveFlag;
parameters.plotResults = plotResults; 
parameters.save_path = save_path; %In the case you'd like to save from within a function
disp('Data Imported')

% Test parameters
parameters.nTrials = 1; %Number of test initializations
parameters.nNets = 1; %Number of network initializations
test_n = 10;

% Varied parameters
variedParam(1).name = 'del_G_syn_E_I';
variedParam(1).range =  linspace(1*10^(-10), 20*10^(-9), test_n);
variedParam(2).name = 'n';
variedParam(2).range =  linspace(100, 500, test_n-1);
variedParam(3).name = 'clusters';
variedParam(3).range =  linspace(1, 10, test_n);
variedParam(4).name = 'mnc';
variedParam(4).range =  linspace(1, 10, test_n-1);

% Combine into one parameter vector to pass to parfor function
parameterSets_vec = combvec(variedParam(:).range);

% Exclude cases where mnc>clusters
for i = 1:length(variedParam)
    for j = 1:length(variedParam)
        if isequal(variedParam(i).name, 'mnc') && isequal(variedParam(j).name, 'clusters')
            parameterSets_vec = parameterSets_vec(:,~[parameterSets_vec(i,:)>parameterSets_vec(j,:)]);
        end
    end
end
   
% Save cases and updated parameters
if saveFlag
    save(strcat(save_path,'/variedParam.mat'),'variedParam')
end
disp('Varied Parameters Set')

gcp; % starts parallel pool if not already running

% Run parallel criticality tests
disp('Running Criticality Tests')
tic
parfor ithParamSet = 1:size(parameterSets_vec, 2) %Run through all parameter combinations
    %For each combination run parallelize_parameter_tests_2
    parallelize_parameter_tests_hpcc(parameters, parameterSets_vec, ...
        ithParamSet, variedParam);  
end
runTime = toc;
sprintf('Program Runtime (s) = %.2f',runTime)
