function parameters = set_depedent_parameters(parameters)
% Takes in "parameters" structure of initialized independent parameters and
% sets (or updates) all dependent parameters.

%Interaction constants
parameters.t_steps = parameters.t_max/parameters.dt; %number of timesteps in simulation
parameters.syn_E = parameters.V_syn_E; %vector of the synaptic reversal potential for excitatory connections
parameters.syn_I = parameters.V_syn_I; %vector of the synaptic reversal potential for inhibitory connections
parameters.IES = ceil(parameters.IEI/parameters.dt); %inter-event-steps = the number of steps to elapse between spikes

%Calculate connection probabilites
% parameters.cluster_n = min(parameters.n*2/parameters.clusters,parameters.n); %number of neurons in a cluster (for small n round(n/3), for large n round(n/5)) 
parameters.cluster_n = round((parameters.mnc*parameters.n) / parameters.clusters) ; % Rather than above method, explicitly declare mnc as a parameter
parameters.npairs = parameters.n*(parameters.n-1); %total number of possible neuron connections
parameters.nclusterpairs = parameters.cluster_n*(parameters.cluster_n - 1)*parameters.clusters; %total number of possible intra-cluster connections
parameters.cluster_prob = min(parameters.conn_prob*parameters.npairs/parameters.nclusterpairs,1); %0.2041; %intra-cluster connection probability
parameters.n_I = round((1-parameters.p_E)*parameters.n); %number of inhibitory neurons

end