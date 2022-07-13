function parameters = set_depedent_parameters(parameters)
% Takes in "parameters" structure of initialized independent parameters and
% sets (or updates) all dependent parameters.

%Interaction constants
parameters.t_steps = parameters.t_max/parameters.dt; %number of timesteps in simulation
parameters.syn_E = parameters.V_syn_E; %vector of the synaptic reversal potential for excitatory connections
parameters.syn_I = parameters.V_syn_I; %vector of the synaptic reversal potential for inhibitory connections

if isfield(parameters, 'IEI')
    parameters.IES = ceil(parameters.IEI/parameters.dt); %inter-event-steps = the number of steps to elapse between spikes
end

% To make Wei and Wie equal, set one to nan before calling set_dependent_parameters
if isnan(parameters.del_G_syn_I_E) & ~isnan(parameters.del_G_syn_E_I)
    parameters.del_G_syn_I_E = parameters.del_G_syn_E_I;
elseif isnan(parameters.del_G_syn_E_I) & ~isnan(parameters.del_G_syn_I_E)
    parameters.del_G_syn_E_I = parameters.del_G_syn_I_E;
elseif isnan(parameters.del_G_syn_E_I) & isnan(parameters.del_G_syn_I_E)
    error('Wei and Wie are both nan')
end

%Calculate connection probabilites
parameters.cluster_n = min(round((parameters.mnc*parameters.n) / parameters.clusters),parameters.n) ; % Rather than above method, explicitly declare mnc as a parameter
parameters.npairs = parameters.n*(parameters.n-1); %total number of possible neuron connections
parameters.nclusterpairs = parameters.cluster_n*(parameters.cluster_n - 1)*parameters.clusters; %total number of possible intra-cluster connections
parameters.cluster_prob = min(parameters.conn_prob*parameters.npairs/parameters.nclusterpairs,1); %0.2041; %intra-cluster connection probability
parameters.n_I = round((1-parameters.p_E)*parameters.n); %number of inhibitory neurons
parameters.n_E = parameters.n-parameters.n_I; %number of excitatory neurons

%The following is a test run of a dependent parameter
parameters.del_G_sra = 60*parameters.del_G_syn_E_E - 2*10^(-7);

end