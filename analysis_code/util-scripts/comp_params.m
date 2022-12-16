
%% Add data folder source

if ispc
    addpath('C:\Users\Jordan\Box\Data\RandNet-Data\RandNet param sweeps')
elseif ismac
    addpath(['/Users/jordan/Library/CloudStorage/Box-Box/Data/RandNet-Data/RandNet param sweeps'])
else
    disp('error')
end


%% Compare parameters of a parameter sweep against a base parameter set
% 2022-07-22T18-29, 2022-08-06T22-55, 2022-08-10T21-41, 2022-10-12T09-22, 2022-10-18T02-35


sim1Name = ['results_', '2022-10-18T02-35', '.mat']
sim1 = load(sim1Name, 'parameters', 'pfsim');

sim2Name = ['results_', '2022-07-22T18-29', '.mat']
sim2 = load(sim2Name, 'parameters', 'pfsim');

if 0
    % exampleDir = '/Users/jordan/Documents/GitHub/RandNet/results/randnet_PF_Envs/';
    % exampleDir = 'C:\Users\Jordan\Documents\GitHub\RandNet\results\randnet_PF\'
    %  exampleDir = 'C:\Users\Jordan\Documents\GitHub\RandNet\results\randnet_PF_param_tests\'
    exampleDir = 'C:\Users\Jordan\Documents\GitHub\RandNet\results\randnet_PF\withCorrelations\'
    
    tmp = load([exampleDir, 'parameters.mat']); sim2.parameters = tmp.parameters;
    tmp = load([exampleDir, 'pfsim.mat']); sim2.pfsim = tmp.pfsim;
end


[common, d1, d2] = comp_struct(sim1.parameters, sim2.parameters)
[common, d1, d2] = comp_struct(sim1.pfsim, sim2.pfsim)
