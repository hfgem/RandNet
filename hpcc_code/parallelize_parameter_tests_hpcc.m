function parallelize_parameter_tests_hpcc(parameters, parameterSets_vec, ...
    ithParamSet, variedParam)
    %_________
    %ABOUT: This function runs through a series of commands to test the
    %outputs of a particular parameter set in comparison to a strict set of
    %criteria. This function is to be used in conjunction with
    %randnet_calculator_memOpt.m or randnet_calculator.m
    %INPUTS:
    %   parameters = a structure that contains the following (only
    %   relevant listed below):
    %       n = Number of neurons in the network
    %       clusters = Number of clusters of neurons in network
    %       mnc = Mean number of clusters to which neurons belong
    %       t_max = maximum time of simulation (s)
    %       dt = timestep of simulation (s)
    %       saveFlag = flag of whether to save results (1 to save, 
    %                   0 otherwise)
    %       save_path = where to save results if desired
    %       plotResults = flag of whether to plot results (1 to plot, 
    %                   0 otherwise)
    %       p_I = probability of an I cell connecting to any other cell
    %       nNets = number of network initializations to test
    %       nTrials = number of simulation initializations to test   
    %       E_events_only = flag to analyze only excitatory neuron behavior
    %       check_criticality = flag to check chaos in network. If == 0, then
    %                   won't check, and the 5th value of the avg_mat will 
    %                   be 0 by default.
    %   parameterSets_vec = a matrix of that contains parameter values to
    %       test: rows = number of parameters, columns = number of sets.
    %   ithParamSet = index of which parameter set to test
    %   variedParam = structure containing the names of parameters being
    %       modified
    %OUTPUTS:
    %   avg_mat = A vector representing average values from all network
    %   initializations and the firing initializations for each network
    %   structure. Only successful tests are averaged:
    %       1. average fraction of neurons spiking in an event
    %       2. average firing rate
    %       3. average event length
    %       4. average number of identified events
    %       5. critical (=1) or not(=0)
    %   allResults = Cell structure storing all results
    %_________
    
%     % For debugging only, output which parameter set is in progress:
%     disp(['Parameter set ', num2str(ithParamSet), '/', num2str(size(parameterSets_vec, 2)), ' started'])

    % Set up parameter values for current parameter set
    for i = 1:size(variedParam, 2)
        parameters.(variedParam(i).name) = parameterSets_vec(i,ithParamSet);
    end
    
    % Update any parameters that are dependent on a varied parameter
    parameters = set_depedent_parameters(parameters);

    %Run network initialization code
    seed = 1; %Random number generator seed for initialization
    network = create_clusters(parameters, 'seed', seed, 'include_all', parameters.include_all, 'global_inhib', parameters.global_inhib);
    %Create input conductance variable
    G_in = single(zeros(parameters.n, parameters.t_steps+1));
    for k = 2:(parameters.t_steps+1)
        G_in(:,k) = G_in(:,k-1)*exp(-parameters.dt/parameters.tau_syn_E);
        G_in(:,k) = G_in(:,k) + parameters.W_gin * [randi(10000, parameters.n, 1)/10000 < (parameters.dt*parameters.rG)];
    end
    clear k
    parameters.('G_in') = G_in;
    clear G_in
    %Create Storage Variables Based on Calculator Code Used
    V_m = parameters.V_reset + randn([parameters.n,1])*parameters.V_m_noise; %set all neurons to baseline reset membrane potential with added noise
    spikes_V_m = randnet_calculator_memOpt(parameters, seed, network, V_m);
    parameters = rmfield(parameters, 'G_in');
    clear seed V_m network
    % count avalanches
    [av_lengths, av_counts] = pull_avalanches(spikes_V_m);
    %Test network for chaotic activity
    test_criticality_hpcc(av_lengths, av_counts, parameters, ithParamSet);
    %Read out time of parameter set completion
    cur_time = clock;
    h_t = string(cur_time(4));
    m_t = string(cur_time(5));
    if length(m_t) == 1
        m_t = '0' + m_t;
    end
    disp(strcat('Time = ',h_t,':',m_t,': Parameter set ', num2str(ithParamSet), '/', num2str(size(parameterSets_vec, 2)), ' complete'))
end % Function