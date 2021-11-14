% Testing the effect of one neuron on another and finding the right
% parameter space to match it

%Load network structure and parameters
load_path = uigetdir('/Users/hannahgermaine/Documents/PhD/','Select Network Load Folder'); %Have user input where they're pulling parameters from
load(strcat(load_path,'/network.mat'))
slashes = find(load_path == '/');
param_path = load_path(1:slashes(end));
load(strcat(param_path,'/parameters.mat'))

%% Test deflections

parameters_copy = parameters;

runs = 20; %number of times to test 

%Change simulation time and initialization type
parameters_copy.t_max = 0.1;
parameters_copy.type = 'neuron';
parameters_copy.connectivity_gain = 0;
parameters_copy.I_coeff = 0;

%Storage
V_m_tests_pre = zeros(runs,parameters_copy.t_max/parameters_copy.dt + 1);
V_m_tests_post = zeros(runs,parameters_copy.t_max/parameters_copy.dt + 1);

for i = 1:runs
    %Create a network copy and keep only one pair of connections
    rng(i)
    network_copy = network;
    conns_ind = find(network_copy.conns);
    conns_ind = conns_ind(randperm(length(conns_ind)));
    network_copy.conns(conns_ind(2:end)) = 0;
    [pre_i,post_i] = find(network_copy.conns);
    %____________________________________
    %___Calculate Dependent Parameters___
    %____________________________________
    cluster_n = min(parameters_copy.n*2/parameters_copy.clusters,parameters.n); %number of neurons in a cluster (for small n round(n/3), for large n round(n/5)) 
    parameters_copy.('cluster_n') = cluster_n;

    %Interaction constants
    t_steps = parameters_copy.t_max/parameters_copy.dt; %number of timesteps in simulation
    E_syn_E = parameters_copy.V_syn_E*ones(parameters_copy.n,1); %vector of the synaptic reversal potential for excitatory connections
    E_syn_I = parameters_copy.V_syn_I*ones(parameters_copy.n,1); %vector of the synaptic reversal potential for inhibitory connections
    IES = ceil(parameters_copy.IEI/parameters_copy.dt); %inter-event-steps = the number of steps to elapse between spikes
    %save for easy calculations
    parameters_copy.('t_steps') = t_steps;
    parameters_copy.('E_syn_E') = E_syn_E;
    parameters_copy.('E_syn_I') = E_syn_I;
    parameters_copy.('IES') = IES;

    %Adding an input current to all cells (one of the options must be uncommented)
    x_in = [0:parameters_copy.dt:parameters_copy.t_max];
    % %Rhythmic current input: (uncomment if desired)
    % I_coeff = 0; %5.1*10^(-10); %set to 0 for no input current
    % I_Hz = 1; %frequency of input - 1 Hz or 60 bpm, 4-8 Hz for theta rhythm
    % I_in = I_coeff*(0.5 + 0.5*sin(I_Hz*2*pi*x_in)); %Approximately theta wave input current
    %Noisy input: (uncomment if desired)
    I_in = parameters_copy.I_coeff*randn(parameters_copy.n,parameters_copy.t_steps+1)*parameters_copy.I_scale; %Generally use -0.5-0.5 nA stimulus
    %save for easy calculations
    parameters_copy.('x_in') = x_in;
    parameters_copy.('I_in') = I_in;

    %Calculate connection probabilites
    npairs = parameters_copy.n*(parameters_copy.n-1); %total number of possible neuron connections
    nclusterpairs = parameters_copy.cluster_n*(parameters_copy.cluster_n - 1)*parameters_copy.clusters; %total number of possible intra-cluster connections
    cluster_prob = min(parameters_copy.conn_prob*npairs/nclusterpairs,1); %0.2041; %intra-cluster connection probability
    p_I = 1 - parameters_copy.p_E; %probability of an inhibitory neuron
    n_I = round(p_I*parameters_copy.n); %number of inhibitory neurons
    %save for easy calculations
    parameters_copy.('npairs') = npairs;
    parameters_copy.('nclusterpairs') = nclusterpairs;
    parameters_copy.('cluster_prob') = cluster_prob;
    parameters_copy.('p_I') = p_I;
    parameters_copy.('n_I') = n_I;
    
    %Create Storage Variables
    I_syn = zeros(parameters_copy.n,parameters_copy.t_steps+1); %synaptic current emitted by each neuron at each timestep (A)
    %synaptic conductances for each neuron at each timestep
    G_syn_I = zeros(parameters_copy.n,parameters_copy.t_steps+1); %conductance for presynaptic inhibitory (S)
    G_syn_E = zeros(parameters_copy.n,parameters_copy.t_steps+1); %conductance for presynaptic excitatory (S)
    V_m = zeros(parameters_copy.n,parameters_copy.t_steps+1); %membrane potential for each neuron at each timestep
    V_m(:,1) = parameters_copy.V_reset + randn([parameters_copy.n,1])*(10^(-3)); %set all neurons to baseline reset membrane potential with added noise
    G_sra = zeros(parameters_copy.n,parameters_copy.t_steps+1); %refractory conductance for each neuron at each timestep
    
    %Now run simulation
    [V_m, ~, ~, ~, ~] = lif_sra_calculator_postrotation(...
            parameters_copy, 1, network_copy, I_syn, G_syn_I, G_syn_E, V_m, G_sra);
    
    V_m_tests_post(i,:) = V_m(post_i,:);
    V_m_tests_pre(i,:) = V_m(pre_i,:);
end

%% Calculate values of the simulated pairs

V_m_max = max(V_m_tests_post,[],2); %maximum V_m value in V
V_m_max_delta = abs(V_m_max - V_m_tests_post(:,1)); %maximum deflection in V
keeps = V_m_max_delta >= 0.001;
V_m_max_avg = mean(V_m_max_delta(keeps)); %average deflection in V
[~,V_m_max_delay] = find(V_m_tests_post == V_m_max); %delay to max in seconds
V_m_delay_avg = mean(V_m_max_delay(keeps))*parameters_copy.dt; %average delay to max in seconds

%V_m_delay_avg = 0.0231
%V_m_max_avg = 0.0061

%% Now test the parameter space of different capacitances and conductances
%to see which match the original values

C_m_vals = linspace(0.1*10^(-12),1*10^(-9),10);
del_G_vals = linspace(0.5*10^(-9),8*10^(-9),10);

parameters = parameters_copy;
parameters_copy = parameters;

V_m_delay_avg_tests = zeros(10,10);
V_m_max_avg_tests = zeros(10,10);

for C_m = 1:10
    %Update params
    parameters_copy.C_m = C_m_vals(C_m);
    for del_G = 1:10
        %Update Params
        parameters_copy.del_G_syn_E = del_G_vals(del_G);
        parameters_copy.del_G_syn_I = del_G_vals(del_G);
        
        %Storage
        V_m_tests_pre = zeros(runs,parameters_copy.t_max/parameters_copy.dt + 1);
        V_m_tests_post = zeros(runs,parameters_copy.t_max/parameters_copy.dt + 1);

        %Run multiple inits
        for i = 1:runs
            %Select neuron pair
            %Create a network copy and keep only one pair of connections
            rng(i)
            network_copy = network;
            conns_ind = find(network_copy.conns);
            conns_ind = conns_ind(randperm(length(conns_ind)));
            network_copy.conns(conns_ind(2:end)) = 0;
            [pre_i,post_i] = find(network_copy.conns);
            
            %Create Storage Variables
            I_syn = zeros(parameters_copy.n,parameters_copy.t_steps+1); %synaptic current emitted by each neuron at each timestep (A)
            %synaptic conductances for each neuron at each timestep
            G_syn_I = zeros(parameters_copy.n,parameters_copy.t_steps+1); %conductance for presynaptic inhibitory (S)
            G_syn_E = zeros(parameters_copy.n,parameters_copy.t_steps+1); %conductance for presynaptic excitatory (S)
            V_m = zeros(parameters_copy.n,parameters_copy.t_steps+1); %membrane potential for each neuron at each timestep
            V_m(:,1) = parameters_copy.V_reset + randn([parameters_copy.n,1])*(10^(-3)); %set all neurons to baseline reset membrane potential with added noise
            G_sra = zeros(parameters_copy.n,parameters_copy.t_steps+1); %refractory conductance for each neuron at each timestep

            %Now run simulation
            [V_m, ~, ~, ~, ~] = lif_sra_calculator_postrotation(...
                    parameters_copy, 1, network_copy, I_syn, G_syn_I, G_syn_E, V_m, G_sra);
            
            V_m_tests_post(i,:) = V_m(post_i,:);
            V_m_tests_pre(i,:) = V_m(pre_i,:);
            
        end
        
        %Calculate V_m average changes
        V_m_max = max(V_m_tests_post,[],2); %maximum V_m value in V
        V_m_max_delta = V_m_max - V_m_tests_post(:,1); %maximum deflection in V
        V_m_max_avg = mean(V_m_max_delta); %average deflection in V
        [~,V_m_max_delay] = find(V_m_tests_post == V_m_max); %delay to max in seconds
        V_m_delay_avg = mean(V_m_max_delay)*parameters_copy.dt; %average delay to max in seconds

        %Store
        V_m_delay_avg_tests(C_m,del_G) = V_m_delay_avg;
        V_m_max_avg_tests(C_m,del_G) = V_m_max_avg;
        
    end
end    
    
%% Plot and compare

%V_m_delay_avg = 0.0231
%V_m_max_avg = 0.0061

%Plot parameter space
figure;
ax1 = subplot(1,2,1);
scatter3(reshape(C_m_vals'.*ones(10,10),[],1),reshape(del_G_vals.*ones(10,10),[],1),reshape(V_m_delay_avg_tests,[],1))
hold on
scatter3(parameters.C_m,parameters.del_G_syn_E,V_m_delay_avg)
xlabel('C_m')
ylabel('del_G')
zlabel('Max Delay in Seconds')
ax1 = subplot(1,2,2);
scatter3(reshape(C_m_vals'.*ones(10,10),[],1),reshape(del_G_vals.*ones(10,10),[],1),reshape(V_m_max_avg_tests,[],1))
hold on
scatter3(parameters.C_m,parameters.del_G_syn_E,V_m_max_avg)
xlabel('C_m')
ylabel('del_G')
zlabel('Average Change in V_m from Baseline to Max (V)')

% Calculate closest values
dist_delay = abs(V_m_delay_avg_tests - V_m_delay_avg);
dist_max = abs(V_m_max_avg_tests - V_m_max_avg);
figure;
subplot(1,2,1)
imagesc(dist_delay)
title('Distance in Delay Value from Original to New')
ylabel('C_m')
xlabel('del_G')
colorbar()
subplot(1,2,2)
imagesc(dist_max)
title('Distance in \Delta Max Value from Original to New')
ylabel('C_m')
xlabel('del_G')
colorbar()