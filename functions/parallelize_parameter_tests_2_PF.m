function [avg_mat, allResults, PFresults] = parallelize_parameter_tests_2(parameters, pfsim, num_nets,...
    num_inits, parameterSets_vec, ithParamSet, variedParam)
    %_________
    %ABOUT: This function runs through a series of commands to test the
    %outputs of a particular parameter set in comparison to a strict set of
    %criteria. This function is to be used in conjunction with
    %network_tests2.m where a parallelized for loop calls this function,
    %parallelize_networks.m where network initializations are parallelized,
    %and parallelize_network_tests.m where network tests are parallelized.
    %
    %INPUTS:
    %   parameters = a structure that contains the following (only
    %   relevant listed below):
    %       n = Number of neurons in the network
    %       clusters = Number of clusters of neurons in network
    %       t_max = maximum time of simulation (s)
    %       dt = timestep of simulation (s)
    %       tau_syn_E = AMPA synaptic decay time constant (s) [Ignoring NMDA as slow and weak]
    %       tau_syn_I = GABA synaptic decay time constant (s)
    %       tau_stdp = STDP decay time constant (s)
    %       E_K = Potassium reversal potential (V)
    %       E_L = Leak reversal potential (V)
    %       G_L = Leak conductance (S) - 10-30 nS range
    %       C_m = Total membrane capacitance (F)
    %       V_m_noise = Magnitude of membrane potential simulation noise (V)
    %       V_th = The threshold membrane potential (V)
    %       V_reset = The reset membrane potential (V)
    %       V_syn_E = Excitatory synaptic reversal potential (V)
    %       V_syn_I = Inhibitory synaptic reversal potential (V)
    %       del_G_syn_E_E = Synaptic conductance step for 
    %               excitatory-excitatory neuron connections following
    %               spike (S)
    %       del_G_syn_E_I = Synaptic conductance step for 
    %               excitatory-inhibitory neuron connections following
    %               spike (S)
    %       del_G_syn_I_I = Synaptic conductance step for 
    %               inhibitory-inhibitory neuron connections following
    %               spike (S)
    %       del_G_syn_I_E = Synaptic conductance step for 
    %               inhibitory-excitatory neuron connections following
    %               spike (S)
    %       del_G_sra = spike rate adaptation conductance step following spike 
    %               ranges from 1-200 *10^(-9) (S)
    %       tau_sra = Spike rate adaptation time constant (s)
    %       connectivity_gain = Amount to increase or decrease connectivity by 
    %               with each spike (more at the range of 1.002-1.005) -
    %               keep at 1 to ensure no connectivity change
    %       G_coeff = input conductance coefficient (setting strength) (S)
    %       G_scale = input conductance scale (ex. nano = 1*10^(-9)) (S)
    %       t_steps = The number of timesteps in the simulation
    %       syn_E = An [n x 1] vector of the synaptic reversal potential for
    %               excitatory connections (V)
    %       syn_I = An [n x 1] vector of the synaptic reversal potential for
    %               inhibitory connections (V)     
    %   num_nets = number of network structures to test per parameter set 
    %   num_inits = number of network initializations to test per network
    %       structure
    %   parameter_vec = a matrix of size [num_params,test_n] that
    %       contains parameter values to test
    %   test_n = number of parameter values to test
    %   ind = what index of parameter combinations is being tested
    %   save_path = where to save results
    %OUTPUTS:
    %   avg_mat = A vector representing average values from all network
    %   initializations and the firing initializations for each network
    %   structure. Namely, the first dimension length is the number of
    %   networks, the second dimension length is the number of
    %   initializations, and the third dimension length is the 3 parameters
    %   from each combination:
    %       1. number of spiking neurons
    %       2. average firing rate
    %       3. average event length
    %_________
    

    % Set up parameter values for current parameter set
    for i = 1:size(variedParam, 2)
        parameters.(variedParam(i).name) = parameterSets_vec(i,ithParamSet);
        pfsim.(variedParam(i).name) = parameterSets_vec(i,ithParamSet);
    end

    
    % Update any parameters that are dependent on a varied parameter
    parameters = set_depedent_parameters(parameters);
    pfsim = set_depedent_parameters(pfsim);

    %Run network initialization code
    resp_mat = zeros(num_nets, 4);
    allResults = cell(1, num_nets) ;
    PFresults = cell(1, num_nets) ;
    for ithNet = 1:num_nets
        
        network = create_clusters(parameters, 'seed', ithNet, 'include_all', parameters.include_all, 'global_inhib', parameters.global_inhib);
        
        mat = zeros(num_inits,4);
        network_spike_sequences = struct; 
        
        %% Run PF simulation
        G_in_PFs = zeros(parameters.n, numel(pfsim.t), pfsim.nEnvironments, pfsim.nTrials);
        for ithEnv = 1:pfsim.nEnvironments
            
            PFresults{ithNet}{ithEnv}.E_indices = network.E_indices;

                    
            for ithTrial = 1:pfsim.nTrials
                %rng(ithTrial)
                % G_in_PFs(:,1,ithEnv,ithTrial) = 1/10* dI(:,ithEnv) * 2*parameters.rGmax * parameters.tau_syn_E + sqrt(1/2*parameters.tau_syn_E*dI(:,ithEnv).^2*2*parameters.rGmax).*randn(parameters.n, 1) ; 
                G_in_PFs(:,1,ithEnv,ithTrial) = zeros(parameters.n, 1) ; 
                for i = 2:numel(pfsim.t)
                        % Exponential decay from last time step
                        G_in_PFs(:,i,ithEnv,ithTrial) = G_in_PFs(:,i-1,ithEnv,ithTrial)*exp(-parameters.dt/parameters.tau_syn_E);
                        % Add spikes from each input source
                        G_in_PFs(:,i,ithEnv,ithTrial) = G_in_PFs(:,i,ithEnv,ithTrial) + ...
                                network.spatialInput{1} .* (1-parameters.PFcontextFrac).*[rand(parameters.n, 1) < (parameters.dt* (parameters.rG * (i/numel(pfsim.t)) ))] + ...
                                network.spatialInput{2} .* (1-parameters.PFcontextFrac).* [rand(parameters.n, 1) < (parameters.dt* ( parameters.rG * ((numel(pfsim.t)-i)/numel(pfsim.t)) ) )] + ...
                                network.contextInput(:,ithEnv) .* [parameters.PFcontextFrac.*ismember(network.all_indices, network.E_indices)]' .* [rand(parameters.n, 1) < (parameters.dt*parameters.rG)] + ...
                                network.contextInput(:,ithEnv) .* [parameters.IcueScale_PF.*ismember(network.all_indices, network.I_indices)]' .* [rand(parameters.n, 1) < (parameters.dt*parameters.rG)]   ;
                end
            end

            opV = zeros(parameters.n, numel(pfsim.t), pfsim.nEnvironments, pfsim.nTrials); % Voltage from all sims
            opS = zeros(parameters.n, numel(pfsim.t), pfsim.nEnvironments, pfsim.nTrials); % Spikes from all sims
            for ithEnv = 1:pfsim.nEnvironments
                for ithTrial = 1:pfsim.nTrials

                    % Set up for simulation
                    V_m = zeros(parameters.n,numel(pfsim.t)); %membrane potential for each neuron at each timestep
                    V_m(:,1) = -60e-3 + 1e-3*randn([parameters.n,1]); %set all neurons to baseline reset membrane potential with added noise
                    pfsim.G_in = G_in_PFs(:,:,ithEnv,ithTrial); 
                    trialSeed = randi(10^6);

                    % PF Simulation
                    [V_m, G_sra, G_syn_E_E, G_syn_I_E, G_syn_E_I, G_syn_I_I, conns] = randnet_calculator(pfsim, trialSeed, network, V_m);
                    opV(:,:,ithEnv,ithTrial) = V_m;
                    opS(:,:,ithEnv,ithTrial) = logical( V_m>parameters.V_th);
                    clear V_m 
                    rmfield(pfsim, 'G_in');
                end
            end
        end
        % keyboard

        % Calculate Place fields and PF-objective score
        allScores = zeros(size(pfsim.nEnvironments));
        
        for ithEnv = 1:pfsim.nEnvironments
            % [linfields, PFpeaksSequence, PFmat] = calculate_linfields(opS, pfsim, pfsim, network, true);
            [~, ~, PFmat] = calculate_linfields(opS, pfsim, pfsim, network, true);
            PFresults{ithNet}{ithEnv}.linfields = PFmat;

            %{
            figure; plot(pfsim.t, V_m(network.E_indices(1:3),:)); ylabel('Vm (V)'); xlabel('Time (s)'); 
            figure; plot(pfsim.t, V_m(network.I_indices(1:3),:)); ylabel('Vm (V)'); xlabel('Time (s)'); 
            figure; plot(pfsim.t, G_in_PFs(1:4,:,1,1)); ylabel('G in (S)'); xlabel('Time (s)');      
            figure; plotSpikeRaster( logical( [ opS(network.E_indices,:,1,1), opS(network.I_indices,:,1,1) ]) ...
            'TimePerBin', parameters.dt, 'PlotType', 'scatter'); % 
            ylabel('Cell');
            figure; plotSpikeRaster( logical(opS(network.I_indices,:,1,1)), ...
            'TimePerBin', parameters.dt, 'PlotType', 'scatter'); % 
            ylabel('Cell');
            %}
            if pfsim.PFscoreFlag % Calculating the score is computationally expensive (due to fitting gaussian curves)
                allScores(i) = calculate_linfieldsScore(linfields, pfsim, pfsim, network)
                %disp(['Env: ', num2str(i), ', Score: ', num2str(allScores(i))])
            end
        end
        if pfsim.PFscoreFlag
            score = mean(allScores);
        end
        %disp(['Mean score: ', num2str(score)])

        %figure; plotSpikeRaster( logical( [ opS(network.E_indices,:,1,1); opS(network.I_indices,:,1,1) ]), 'TimePerBin', parameters.dt, 'PlotType', 'scatter');
        %figure; histogram( sum(opS, [2:4])./parameters.t_max./parameters.nTrials, 50 )
        
        clear G_in_PFs opV opS
        
        %% Run preplay simulation
        for ithTest = 1:num_inits
            seed = ithTest;

            %Create input conductance variable
            if parameters.usePoisson
                G_in = single(zeros(parameters.n, parameters.t_steps+1));
                
                for k = 2:(parameters.t_steps+1)
                    G_in(:,k) = G_in(:,k-1)*exp(-parameters.dt/parameters.tau_syn_E);

                    % G_in(:,k) = G_in(:,k) + network.contextInput .* [rand(parameters.n, 1) < (parameters.dt*parameters.rG)];
                    G_in(:,k) = G_in(:,k) + [network.contextInput .* [1     .*              ismember(network.all_indices, network.E_indices)]' .* [rand(parameters.n, 1) < (parameters.dt*parameters.rG)] + ...
                                             network.contextInput .* [parameters.IcueScale.*ismember(network.all_indices, network.I_indices)]'  .* [rand(parameters.n, 1) < (parameters.dt*parameters.rG)]] ;
                end
            else
                G_in = (parameters.G_std*randn(parameters.n,parameters.t_steps+1))+parameters.G_mean;
                G_in(G_in<0) = 0;
            end
            parameters.('G_in') = G_in;
            clear G_in
            
            %Run model
            
            %Create Storage Variables
            % V_m = zeros(parameters.n,parameters.t_steps+1); %membrane potential for each neuron at each timestep
            % V_m(:,1) = parameters.V_reset + randn([parameters.n,1])*parameters.V_m_noise; %set all neurons to baseline reset membrane potential with added noise
            % [V_m, ~, ~, ~, ~] = randnet_calculator(parameters, seed, network, V_m);
            % E_spikes_V_m = sparse(V_m(network.E_indices,:) >= parameters.V_th);
            V_m = parameters.V_reset + randn([parameters.n,1])*parameters.V_m_noise; %set all neurons to baseline reset membrane potential with added noise
            spikeMat = randnet_calculator_memOpt(parameters, seed, network, V_m);
            parameters = rmfield(parameters, 'G_in');
            E_spikes_V_m = spikeMat(network.E_indices,:); clear spikeMat
            

            % detect events and compute outputs
            % [network_spike_sequences, network_cluster_sequences, outputVec] = detect_events(parameters, network, V_m , j, network_spike_sequences, network_cluster_sequences);

            [trialResults] = detect_PBE(E_spikes_V_m, parameters);
            if ithTest == 1 % append trialResults struct to network results struct
                network_spike_sequences = trialResults;
            else
                network_spike_sequences = [network_spike_sequences, trialResults]; 
            end
            
            % Calculate mean n spikes within events
            overallmeanspikes = 0;
            for ithEvent = 1:size(trialResults.events, 1)
                eventSpikes = E_spikes_V_m(:,trialResults.events(ithEvent,1):trialResults.events(ithEvent,2));
                nSpikes = sum(eventSpikes, 2);
                meannSpikes = sum(nSpikes)/sum(nSpikes>0);
                overallmeanspikes = overallmeanspikes+ meannSpikes;
            end
            overallmeanspikes = overallmeanspikes./size(trialResults.events, 1);
            
            % Overall simulation statistics
            allResults{ithNet}{ithTest}.ithInit = ithTest;
            allResults{ithNet}{ithTest}.numEvents = numel(network_spike_sequences(ithTest).event_lengths); % number of detected events
            allResults{ithNet}{ithTest}.fracFire =  mean(sum(E_spikes_V_m, 2)>0); % Fraction of cells that fire at all during simulation
            allResults{ithNet}{ithTest}.frac_participation = mean([network_spike_sequences(ithTest).frac_spike{:}]); % mean fraction of cells firing per event
            allResults{ithNet}{ithTest}.meanRate = mean(sum(E_spikes_V_m, 2)/parameters.t_max); % mean over cells' average firing rate
            allResults{ithNet}{ithTest}.stdRate = std(sum(E_spikes_V_m, 2)/parameters.t_max); % STD over cells' average firing rate
            allResults{ithNet}{ithTest}.meanCellnEventSpikes = overallmeanspikes; % mean n event spikes over cells that spike


            % Stats for each detected event
            allResults{ithNet}{ithTest}.eventLength = network_spike_sequences(ithTest).event_lengths; % duration in seconds of all detected events
            allResults{ithNet}{ithTest}.eventParticipation = [network_spike_sequences(ithTest).frac_spike{:}]; % fraction of cells that fired in each event
                        
            % Save ranks_vec
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