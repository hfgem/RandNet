
%% Add data folder source

if ispc
    addpath('C:\Users\Jordan\Box\Data\RandNet-Data\RandNet param sweeps')
elseif ismac
    addpath(['/Users/jordan/Library/CloudStorage/Box-Box/Data/RandNet-Data/RandNet param sweeps'])
else
    disp('error')
end


%% Compare parameters of a parameter sweep against a base parameter set

sourceSweep = 'results_2022-08-06T22-55.mat'
sweepParameters = load(sourceSweep, 'parameters')
sweepPFsim = load(sourceSweep, 'pfsim')

load('/Users/jordan/Documents/GitHub/RandNet/results/randnet_PF_Envs/parameters.mat')
load('/Users/jordan/Documents/GitHub/RandNet/results/randnet_PF_Envs/pfsim.mat')


[common, d1, d2] = comp_struct(sweepParameters.parameters, parameters)
[common, d1, d2] = comp_struct(sweepPFsim.pfsim, pfsim)
