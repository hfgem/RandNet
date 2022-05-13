function [avg_mat, allResults] = parallelize_parameter_tests_2(parameters,...
    parameterSets_vec, ithParamSet, variedParam, optFlag)
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
    %   parameterSets_vec = a matrix of that contains parameter values to
    %       test: rows = number of parameters, columns = number of sets.
    %   ithParamSet = index of which parameter set to test
    %   variedParam = structure containing the names of parameters being
    %       modified
    %   optFlag = flag to use optimized randnet_calculator_memOpt.m code
    %       or unoptimized randnet_calculator.m code.
    %OUTPUTS:
    %   avg_mat = A vector representing average values from all network
    %   initializations and the firing initializations for each network
    %   structure. Only successful tests are averaged:
    %       1. average fraction of neurons spiking in an event
    %       2. average firing rate
    %       3. average event length
    %       4. average number of identified events
    %_________

    % Set up parameter values for current parameter set
    for i = 1:size(variedParam, 2)
        parameters.(variedParam(i).name) = parameterSets_vec(i,ithParamSet);
    end
    
    % Update any parameters that are dependent on a varied parameter
    parameters = set_depedent_parameters(parameters);

    %Run network initialization code
    resp_mat = zeros(parameters.nNets, 4);
    allResults = cell(1, parameters.nNets) ;
    for ithNet = 1:parameters.nNets
        
        network = create_clusters(parameters, 'seed', ithNet, 'include_all', parameters.include_all, 'global_inhib', parameters.global_inhib);
        
        if parameters.saveFlag == 1
            if ~isfolder(strcat(parameters.save_path,'/networks'))
                mkdir(parameters.save_path,'/networks')
            end
            save(strcat(parameters.save_path,'/networks','/network_',string(ithParamSet),'_',string(ithNet),'.mat'),'network')
        end
        
        mat = zeros(parameters.nTrials,4);
        network_spike_sequences = struct; 
        for ithTest = 1:parameters.nTrials
            
            %Random number generator seed for initialization
            seed = ithTest; 
            
            %Create input conductance variable
            if parameters.usePoisson
                G_in = single(zeros(parameters.n, parameters.t_steps+1));
                for k = 2:(parameters.t_steps+1)
                    G_in(:,k) = G_in(:,k-1)*exp(-parameters.dt/parameters.tau_syn_E);
                    G_in(:,k) = G_in(:,k) + parameters.W_gin * [rand(parameters.n, 1) < (parameters.dt*parameters.rG)];
                end
            else
                G_in = (parameters.G_std*randn(parameters.n,parameters.t_steps+1))+parameters.G_mean;
                G_in(G_in<0) = 0;
            end
            parameters.('G_in') = G_in;
            clear G_in
            
            %Run model
            
            %Create Storage Variables Based on Calculator Code Used
            if optFlag == 1 %Use randnet_calculator_memOpt.m
                V_m = parameters.V_reset + randn([parameters.n,1])*parameters.V_m_noise; %set all neurons to baseline reset membrane potential with added noise
                spikes_V_m = randnet_calculator_memOpt(parameters, seed, network, V_m);
                parameters = rmfield(parameters, 'G_in');
            else %Use randnet_calculator.m
                V_m = zeros(parameters.n,parameters.t_steps+1); %membrane potential for each neuron at each timestep
                V_m(:,1) = parameters.V_reset + randn([parameters.n,1])*parameters.V_m_noise; %set all neurons to baseline reset membrane potential with added noise
                [V_m, ~, ~, ~, ~] = randnet_calculator(parameters, seed, network, V_m);
                spikes_V_m = V_m >= parameters.V_th;
            end
            
            %Set which spikes are analyzed based on E_events_only flag
            if parameters.E_events_only == 1 %Use only excitatory neurons in analyses
                used_spikes_mat = spikes_V_m(network.E_indices,:);
            else %Use all neurons in analyses
                used_spikes_mat = spikes_V_m;
            end
            
            % detect events and compute outputs
            if strcmp(parameters.eventType,'PBE')
                [trialResults] = detect_PBE(used_spikes_mat, parameters);
            else
                [trialResults, outputVec] = detect_events(parameters, used_spikes_mat, ithParamSet, ithNet, ithTest);
            end
            
            % append trialResults struct to network results struct
            if ithTest == 1
                network_spike_sequences = trialResults;
            else
                network_spike_sequences = [network_spike_sequences, trialResults];  %#ok<*AGROW>
            end
            
            % Overall simulation statistics
            allResults{ithNet}{ithTest}.ithInit = ithTest;
            allResults{ithNet}{ithTest}.numEvents = numel(network_spike_sequences(ithTest).event_lengths); % number of detected events
            allResults{ithNet}{ithTest}.fracFire =  mean(sum(used_spikes_mat, 2)>0); % Fraction of cells that fire at all during simulation
            allResults{ithNet}{ithTest}.frac_participation = mean([network_spike_sequences(ithTest).frac_spike{:}]); % mean fraction of cells firing per event
            allResults{ithNet}{ithTest}.meanRate = mean(sum(used_spikes_mat, 2)/parameters.t_max); % mean over cells' average firing rate
            allResults{ithNet}{ithTest}.stdRate = std(sum(used_spikes_mat, 2)/parameters.t_max); % STD over cells' average firing rate

            % Stats for each detected event
            allResults{ithNet}{ithTest}.eventLength = network_spike_sequences(ithTest).event_lengths; % duration in seconds of all detected events
            allResults{ithNet}{ithTest}.eventParticipation = [network_spike_sequences(ithTest).frac_spike{:}]; % fraction of cells that fired in each event
                        
            % Save spike_order and ranks_vec
            allResults{ithNet}{ithTest}.spike_order = network_spike_sequences(ithTest).spike_order;
            allResults{ithNet}{ithTest}.ranksVec = network_spike_sequences(ithTest).ranks_vec;
            
            % Main output statistics
            %outputVec(1) = allResults{ithNet}{ithTest}.fracFire; % fraction of spiking neurons over entire simulation
            outputVec(1) = allResults{ithNet}{ithTest}.frac_participation; % fraction of spiking neurons over identified events
            outputVec(2) = allResults{ithNet}{ithTest}.meanRate ; %average firing rate
            outputVec(3) = mean(allResults{ithNet}{ithTest}.eventLength); % Average event length
            outputVec(4) = allResults{ithNet}{ithTest}.numEvents; % Number of events

            mat(ithTest,:) = outputVec; 
            %allResults{ithNet}{j} = allTrialResults;
        end % trial loop
        mat(isnan(mat)) = 0;
        resp_mat(ithNet,:) = sum(mat,1) ./ sum(mat > 0,1); %Only averaging those that did successfully produce data        
        
    end % Network loop
    resp_mat(isnan(resp_mat)) = 0;
    avg_mat = sum(resp_mat,1)./sum(resp_mat > 0,1); %Only averaging those that had results

    disp(['Parameter set ', num2str(ithParamSet), '/', num2str(size(parameterSets_vec, 2)), ' complete'])

end % Function