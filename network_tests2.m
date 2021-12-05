%Testing different network properties for sequence generation

%______ABOUT______:
%This code uses a previously successful set of network parameters by
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
%TRY PINK NOISE - HAVE A DECAY AND TIME CONSTANT THAT KEEPS CONDUCTANCE
%SOMEWHAT 'SMOOTH'
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
test_n = 11;

%Define parameter vectors - and which parameters to test. Do not forget to
%modify parallelize_parameter_tests.m
del_G_sra_vec = linspace(200*10^(-9),400*10^(-9),test_n);
G_coeff_vec = linspace(500,1000,test_n);

parameter_vec = [del_G_sra_vec; G_coeff_vec];
clear del_G_sra_vec G_i_vec

%Save parameter values
save(strcat(save_path,'/parameter_vec.mat'),'parameter_vec','-v7.3')

%Set up storage matrix
%success = zeros(test_n,test_n,test_n,test_n);
[num_params, ~] = size(parameter_vec);
success = zeros(test_n*ones(1,num_params));
%dim 1 = del_G_sra
%dim 2 = G_i

%Test parameters
num_nets = 5; %number of network structure initializations
num_inits = 5; %number of firing initializations

%% Run Grid Search
%Start parallel pool for parallelizing the grid search
% parpool(4)

%Loop through all parameter sets
parfevalOnAll(gcp(), @warning, 0, 'off', 'MATLAB:singularMatrix');
parfor ind = 1:test_n^2, success(ind) = parallelize_parameter_tests(parameters,num_nets,...
    num_inits, parameter_vec, test_n, ind); end
save(strcat(save_path,'/success.mat'),'success','-v7.3')

%% Visualize Grid Search Results
%Recall dimensions:
%dim 1 = tau_sra
%dim 2 = del_G_sra
%dim 3 = del_G_syn_E
%dim 4 = del_G_syn_I

%__________
%Calculate and plot average values as function of parameter values to find
%low dimensional trends
%__________

avg_tau_sra = squeeze(sum(success,[2,3,4])/(test_n^3));
avg_del_G_sra = squeeze(sum(success,[1,3,4])/(test_n^3));
avg_del_G_syn_E = squeeze(sum(success,[1,2,4])/(test_n^3));
avg_del_G_syn_I = squeeze(sum(success,[1,2,3])/(test_n^3));

%Visualize Average Trends
figure;
ax1 = subplot(2,2,1);
plot(parameter_vec(1,:),avg_tau_sra)
title('Avg \tau_{SRA}')
ax2 = subplot(2,2,2);
plot(parameter_vec(2,:),avg_del_G_sra)
title('Avg \delta G_{SRA}')
ax3 = subplot(2,2,3);
plot(parameter_vec(3,:),avg_del_G_syn_E)
title('Avg \delta G_{syn_E}')
ax4 = subplot(2,2,4);
plot(parameter_vec(4,:),avg_del_G_syn_I)
title('Avg \delta G_{syn_I}')

%__________
%Calculate and plot 2D results for parameter pairings
%__________

%2D SRA Visualization
avg_SRA_success = squeeze(sum(success,[3,4])/test_n^2);
figure;
imagesc(avg_SRA_success)
colorbar()
xticks(1:test_n)
xticklabels(parameter_vec(2,:))
yticks(1:test_n)
yticklabels(parameter_vec(1,:))
xlabel('\delta G_{SRA}')
ylabel('\tau_{SRA}')
title('Avg SRA Success')

%2D Synaptic Visualization
avg_syn_success = squeeze(sum(success,[1,2])/test_n^2);
figure;
imagesc(avg_syn_success)
colorbar()
xticks(1:test_n)
xticklabels(parameter_vec(4,:))
yticks(1:test_n)
yticklabels(parameter_vec(3,:))
xlabel('\delta G_{syn_I}')
ylabel('\delta G_{syn_E}')
title('Avg Synaptic Step Success')

%% Custom 3D Visualizations Based on Data
%Recall dimensions:
%dim 1 = tau_sra
%dim 2 = del_G_sra
%dim 3 = del_G_syn_E
%dim 4 = del_G_syn_I

%3D visual of success when del_G_syn_E is set
cmap = jet(100);
colormap(cmap)
x_inds = reshape(repmat([1:test_n]'.*ones(test_n,test_n),1,1,test_n),[],1);
y_inds = reshape(repmat([1:test_n].*ones(test_n,test_n),1,1,test_n),[],1);
z_inds = reshape(permute(repmat([1:test_n]'.*ones(test_n,test_n),1,1,test_n),[3,2,1]),[],1);    
for i = 1:test_n
    figure;
    success_set_syn_E = squeeze(success(:,:,i,:));
    color_mappings = reshape(round(success_set_syn_E*100),[],1);
    color_mappings(color_mappings == 0) = 1;
    cmap_vals = cmap(color_mappings,:);
    scatter3(x_inds,y_inds,z_inds,[],cmap_vals,'filled')
    xticks(test_n)
    xticklabels(parameter_vec(1,:))
    xlabel('\tau_{SRA}')
    yticks(test_n)
    yticklabels(parameter_vec(2,:))
    ylabel('\delta G_{SRA}')
    zticks(test_n)
    zticklabels(parameter_vec(3,:))
    zlabel('\delta G_{syn_I}')
    colormap(cmap)
    colorbar()
end
clear i success_set_syn_E color_mappings cmap_vals x_inds y_inds z_inds

%2D visual of success when del_G_syn_E is set
figure;
sub_sqrt = ceil(sqrt(test_n));
axes = [];
colormap(jet)
for i = 1:test_n
    del_E_syn_val = parameter_vec(3,i);
    del_I_syn_exp = del_E_syn_val*1.4;
    [~, del_I_syn_ind] = min(abs(parameter_vec(4,:) - del_I_syn_exp));
    ax = subplot(sub_sqrt,sub_sqrt,i);
    success_mat = squeeze(success(:,:,i,del_I_syn_ind));
    imagesc(success_mat)
    colorbar()
    xticks(1:test_n)
    xticklabels(parameter_vec(1,:))
    yticks(1:test_n)
    yticklabels(parameter_vec(2,:))
    xlabel('\tau_{SRA}')
    ylabel('\Delta_{G_SRA}')
    title(strcat('\Delta_{G_{syn_E}} = ',string(del_E_syn_val)))
    axes = [axes, ax];
end
linkaxes(axes)

%% Plotting the relationship for 2 parameter test
[test_n,~] = size(success); 

%Finding the maximum values and indices
[max_val, max_ind] = max(success,[],1);
all_max_ind = zeros(test_n, test_n);
for i = 1:test_n
all_max_ind(:,i) = [success(:,i) == max_val(i)]';
end
[x_ind, y_ind] = find(all_max_ind);
y_val = parameter_vec(1,x_ind);
x_val = parameter_vec(2,y_ind);
figure;
scatter(x_val, y_val, 'filled')
xlabel('\Delta G_{SRA}')
ylabel('\tau_{SRA}')

%Fit a curve - chose exponential due to apparent plateau and quadratic due
%to curvature.
fit_curve_exp = fit(x_val', y_val', 'exp1');
fit_curve_poly = fit(x_val', y_val', 'poly1');
hold on
plot(unique(x_val), fit_curve_exp(unique(x_val)), 'DisplayName', 'Exponential Fit')
plot(unique(x_val), fit_curve_poly(unique(x_val)), 'DisplayName', 'Polynomial Fit')
legend()
