% STDP Tests 2
% This code allows a user to strengthen a single initialization of a
% sequence through a neuron initialization, then run a long simulation with
% noisy input current to see how many current sequences are highly
% correlated with the "learned" sequence

%% Load Network and Parameters

%Select folder of good network structure and parameters to use in tests
load_path = uigetdir('/Users/hannahgermaine/Documents/PhD/','Select network.m Load Folder'); %Have user input where they'd like network structure to be loaded from

% %_________________________________
% %___Load Independent Parameters___
% %_________________________________
load(strcat(load_path,'/network.mat'))
load(strcat(load_path,'/parameters.mat'))

%% Set Neuron Simulation Parameters

%____Neuron initialization____
neur_seed = 1; %probability of a neuron participating in the initializing group
num_init = 1000; %how many times to run neuron initialization and STDP

%____________________________________
%___Assign STDP Rules Desired________
%____________________________________
%Here you can modify the STDP rules previously saved
parameters.type = 'neuron'; %Each sequence initialization will be through setting neurons to threshold
parameters.tau_stdp = 1*10^(-3); %STDP time constant (s)
parameters.connectivity_gain = 0.005; %amount to increase or decrease connectivity by with each STDP rule (more at the range of 0.002-0.005)

%____________________________________
%___Calculate Dependent Parameters___
%____________________________________
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

%save number of initializations
parameters.('test_val_max') = num_init;

%Adding an input conductance to all cells (one of the options must be uncommented)
parameters.G_coeff = 0;

%% Run Neuron Simulation with STDP

%If this section has already been run, reload a new network.mat
load(strcat(load_path,'/network.mat'))

spikes_struct = struct;
conns_struct = struct;

for i = 1:parameters.test_val_max %Run repeats of initialization
    
    %Create Storage Variables
    %synaptic conductances for each neuron at each timestep
    G_syn_I = zeros(parameters.n,parameters.t_steps+1); %conductance for presynaptic inhibitory (S)
    G_syn_E = zeros(parameters.n,parameters.t_steps+1); %conductance for presynaptic excitatory (S)
    V_m = zeros(parameters.n,parameters.t_steps+1); %membrane potential for each neuron at each timestep
    V_m(:,1) = parameters.V_reset + randn([parameters.n,1])*(10^(-3)); %set all neurons to baseline reset membrane potential with added noise
    G_sra = zeros(parameters.n,parameters.t_steps+1); %refractory conductance for each neuron at each timestep
    
    %Run model
    [V_m, ~, ~, ~, conns] = lif_sra_calculator_postrotation(...
    parameters, neur_seed, network, G_syn_I, G_syn_E, V_m, G_sra);

    %Update connection matrix
    network.conns = conns;
    
    %Find spike profile
    spikes_V_m = V_m >= parameters.V_th;
    [spikes_x,spikes_t] = find(spikes_V_m);
    max_time = max(spikes_t);
    
    %Save results
    if i == 1
        spikes_struct(1).repeat_number = i;
        spikes_struct(1).spikes_Vm = spikes_V_m;
        spikes_struct(1).spikes_x = spikes_x;
        spikes_struct(1).spikes_t = spikes_t;
        spikes_struct(1).max_dt = max_time;
        spikes_struct(1).max_time = max_time*parameters.dt;
        conns_struct(1).repeat_number = i;
        conns_struct(1).conns = conns;
    elseif mod(i,10) == 0
        spikes_struct(end+1).repeat_number = i;
        spikes_struct(end).spikes_Vm = spikes_V_m;
        spikes_struct(end).spikes_x = spikes_x;
        spikes_struct(end).spikes_t = spikes_t;
        spikes_struct(end).max_dt = max_time;
        spikes_struct(end).max_time = max_time*parameters.dt;
        conns_struct(end+1).repeat_number = i;
        conns_struct(end).conns = conns;
    end

end    
clear i G_syn_I G_syn_E V_m G_sra conns spikes_V_m spikes_x spikes_t max_time 
save(strcat(load_path,'/conns_struct.mat'),'conns_struct','-v7.3')
save(strcat(load_path,'/spikes_struct.mat'),'spikes_struct','-v7.3')

%% Visualize connectivity matrix changes

%Visualize first and last sequence
spikes_1 = unique(spikes_struct(1).spikes_x,'stable');
non_spiking_1 = setdiff(1:parameters.n,spikes_1);
spikes_1 = [spikes_1;non_spiking_1'];
spikes_end = unique(spikes_struct(end).spikes_x,'stable');
non_spiking_end = setdiff(1:parameters.n,spikes_end);
spikes_end = [spikes_end;non_spiking_end'];
dt = spikes_struct(1).max_time/spikes_struct(1).max_dt;
max_1 = spikes_struct(1).max_dt;
max_2 = spikes_struct(end).max_dt;
max_time = max(max_1,max_2);
figure;
ax1 = subplot(2,2,1);
imagesc(spikes_struct(1).spikes_Vm(spikes_1,1:max_time))
ticks = xticks();
x_ticklabels = cellstr(string(round(ticks*dt,3)));
xticklabels(x_ticklabels)
title(sprintf('Sequence %i: Ordered by Sequence %i',spikes_struct(1).repeat_number,spikes_struct(1).repeat_number))
xlabel('Time (s)')
ylabel('Neuron Reordered Index')
ax2 = subplot(2,2,2);
imagesc(spikes_struct(1).spikes_Vm(spikes_end,1:max_time))
ticks = xticks();
x_ticklabels = cellstr(string(round(ticks*dt,3)));
xticklabels(x_ticklabels)
title(sprintf('Sequence %i: Ordered by Sequence %i',spikes_struct(1).repeat_number,spikes_struct(end).repeat_number))
xlabel('Time (s)')
ylabel('Neuron Reordered Index')
ax3 = subplot(2,2,3);
imagesc(spikes_struct(end).spikes_Vm(spikes_1,1:max_time))
ticks = xticks();
x_ticklabels = cellstr(string(round(ticks*dt,3)));
xticklabels(x_ticklabels)
title(sprintf('Sequence %i: Ordered by Sequence %i',spikes_struct(end).repeat_number,spikes_struct(1).repeat_number))
xlabel('Time (s)')
ylabel('Neuron Reordered Index')
ax4 = subplot(2,2,4);
imagesc(spikes_struct(end).spikes_Vm(spikes_end,1:max_time))
ticks = xticks();
x_ticklabels = cellstr(string(round(ticks*dt,3)));
xticklabels(x_ticklabels)
title(sprintf('Sequence %i: Ordered by Sequence %i',spikes_struct(end).repeat_number,spikes_struct(end).repeat_number))
xlabel('Time (s)')
ylabel('Neuron Reordered Index')
linkaxes([ax1,ax2,ax3,ax4])
clear ax1 ax2 ax3 ax4 spikes_1 non_spiking_1 spikes_end non_spiking_end max_1 max_2 max_time

%Visualize connectivity matrix changes
spikes_1 = spikes_struct(1).spikes_x;
not_spiked = setdiff(1:parameters.n,spikes_1);
spikes_1 = [spikes_1;not_spiked'];
[~,num_tests_saved] = size(spikes_struct);
inds = round(linspace(1,num_tests_saved,9));
figure;
for i = 1:9
    subplot(3,3,i)
    imagesc(conns_struct(inds(i)).conns(spikes_1,spikes_1))
    colorbar()
    title(sprintf('Repeat = %i',conns_struct(inds(i)).repeat_number))
    xlabel('Index')
    ylabel('Index')
end
sgtitle('Connectivity Matrix')
clear spikes_1 not_spiked num_tests_saved inds i

%% Set Current Initialization Parameters

parameters.type = 'current';
parameters.connectivity_gain = 0; %Turn off STDP

%Input conductance parameters
G_coeff = 16;
parameters.G_coeff = G_coeff;
G_scale = 1*10^(-9);
parameters.G_scale = G_scale;

%Simulation time parameters
t_max = 600; %Time in seconds
parameters.t_max = t_max;
test_val_max = 1;
parameters.test_val_max = test_val_max;
t_steps = parameters.t_max/parameters.dt; %number of timesteps in simulation
parameters.t_steps = t_steps;

%% Run simulation of current initialization

load(strcat(load_path,'/conns_struct.mat'))
curr_seed = 1; %Seed for current initialization
conns_to_test = [1,round(length(conns_struct)/2),length(conns_struct)];
conns_to_test_vals = [conns_struct(conns_to_test).repeat_number];

curr_struct = struct; %to save current sequences
for i = 1:parameters.test_val_max %Run repeats of initialization
    
    %Create Storage Variables
    %synaptic conductances for each neuron at each timestep
    G_syn_I = zeros(parameters.n,parameters.t_steps+1); %conductance for presynaptic inhibitory (S)
    G_syn_E = zeros(parameters.n,parameters.t_steps+1); %conductance for presynaptic excitatory (S)
    V_m = zeros(parameters.n,parameters.t_steps+1); %membrane potential for each neuron at each timestep
    V_m(:,1) = parameters.V_reset + randn([parameters.n,1])*(10^(-3)); %set all neurons to baseline reset membrane potential with added noise
    G_sra = zeros(parameters.n,parameters.t_steps+1); %refractory conductance for each neuron at each timestep
    
    for c = 1:length(conns_to_test)
        %Update connectivity matrix
        network.conns = conns_struct(conns_to_test(c)).conns; %Use the latest connectivity matrix
        
        % Run model
        [V_m, ~, ~, ~, conns] = lif_sra_calculator_postrotation(...
        parameters, curr_seed, network, G_syn_I, G_syn_E, V_m, G_sra);

        % Find spike profile
        spikes_V_m = V_m >= parameters.V_th;
        [spikes_x,spikes_t] = find(spikes_V_m);
        max_time = max(spikes_t);
        spiking_neurons = unique(spikes_x, 'stable');
        spikes_V_m_reordered = spikes_V_m(spiking_neurons,1:max_time);

        % Save basic results
        curr_struct(c).spikes_Vm = spikes_V_m;
        curr_struct(c).spikes_x = spikes_x;
        curr_struct(c).spikes_t = spikes_t;
        curr_struct(c).max_dt = max_time;
        curr_struct(c).max_time = max_time*parameters.dt;
        curr_struct(c).learning_repeats = conns_to_test_vals(c);
        curr_struct(c).conns = conns;
        curr_struct(c).spiking_neurons = spiking_neurons;

        % Run through all events
        if length(spiking_neurons) >= parameters.event_cutoff*parameters.n

            % Find maximum firing rate + average maximum firing rates of neurons
            all_fr = sum(curr_struct(c).spikes_Vm,2)/parameters.t_max; 
            max_fr = max(all_fr);
            avg_fr = mean(all_fr);
            display(avg_fr)

            % TEST 2: The firing rate must fall within a realistic range
            if and(avg_fr>= 0.02, avg_fr <= 1.5)
                % Find event times
                events = []; 
                event_lengths = [];
                last_start = spikes_t(1);
                last_time = spikes_t(1);
                spike_count = 1;
                for t_i = 2:length(spikes_t)
                    s_i = spikes_t(t_i);
                    if s_i - last_time <= parameters.IES
                        last_time = s_i;
                        spike_count = spike_count + 1;
                    else
                        if (last_start ~= last_time) && (spike_count > parameters.event_cutoff*parameters.n) %weed out events w/ too few spikes
                            events = [events; [last_start, last_time]]; %#ok<AGROW> %add the last range of spikes to the events vector
                            event_lengths = [event_lengths, (last_time - last_start)*parameters.dt]; %#ok<*AGROW>
                        end
                        last_start = s_i;
                        last_time = s_i;
                        spike_count = 1;
                    end
                end
                if (last_start ~= last_time) && (spike_count > parameters.event_cutoff*parameters.n) %weed out events w/ too few spikes
                    events = [events; [last_start, last_time]]; %#ok<AGROW> %add the last interval
                    event_lengths = [event_lengths, (last_time - last_start)*parameters.dt]; %#ok<*SAGROW>
                end
                avg_event_length = mean(event_lengths);
                [num_events,~] = size(events);
                curr_struct(c).num_events = num_events;

                for j = 1:num_events
                    curr_struct(c).events(j).times = events(j,:);
                    curr_struct(c).events(j).event_length = event_lengths(j);
                    % TEST 3: The sequence(s) of firing is(are) within
                    % reasonable lengths.
                    if and(avg_event_length >= 0.01, avg_event_length <= 0.15)
                        % store spike orders for each event
                        event_spikes = spikes_V_m(:,events(j,1):events(j,2));
                        [e_spikes_x, ~] = find(event_spikes);
                        spike_order = unique(e_spikes_x,'stable');
                        curr_struct(c).events(j).spike_order = spike_order;
                        % store ranks for each neuron
                        ranks_vec = zeros(1,parameters.n);
                        for k = 1:length(spike_order)
                            n_ind = spike_order(k);
                            ranks_vec(1,n_ind) = k;
                        end
                        curr_struct(c).events(j).spike_ranks = ranks_vec;
                        % store nonspiking neurons
                        nonspiking_neurons = isnan(ranks_vec./ranks_vec);
                        curr_struct(c).events(j).nonspiking_neurons = nonspiking_neurons;
                    end
                    clear e_i event_spikes e_spikes_x spike_order ranks_vec k n_ind nonspiking_neurons
                end %Indiv sequence loop
            end %Firing rate loop
        end %Number of neurons in sequence loop
    end %End of connectivity matrix loop
end
clear events event_lengths last_start last_time spike_count t_i s_i ...
    avg_event_length num_events j i all_fr max_fr avg_fr spikes_x spikes_t ...
    spikes_V_m V_m G_sra G_syn_I G_syn_E conns

save(strcat(load_path,'/curr_struct.mat'),'curr_struct','-v7.3')

%% Correlate the Current Sequences to the Original Neuron Sequence

%load(strcat(load_path,'/curr_struct.mat'))
%load(strcat(load_path,'/spikes_struct.mat'))

[~,num_conns] = size(curr_struct);
[~,num_repeat_tests] = size(spikes_struct);
compare_to = 1; %Index of neuron-initialized learned sequence to compare to

neur_seq = unique(spikes_struct(compare_to).spikes_x,'stable');
non_spiking_neur = setdiff(1:parameters.n,neur_seq);
neur_seq = [neur_seq;non_spiking_neur'];

num_curr_seq = min([curr_struct.num_events]);
corr_vals = zeros(num_conns,num_curr_seq);

for c = 1:num_conns
    for i = 1:num_curr_seq
        curr_seq = unique(curr_struct(c).events(i).spike_order,'stable');
        non_spiking_curr = setdiff(1:parameters.n,curr_seq);
        curr_seq = [curr_seq;non_spiking_curr'];
        corr_vals(c,i) = corr(neur_seq,curr_seq);
    end    
end
clear c i curr_seq non_spiking_curr

colors = ['r','g','b'];

figure;
hold on
for c = 1:num_conns
    histogram(corr_vals(c,:),10,'FaceColor',colors(c),'DisplayName',sprintf('After %i Repeats',curr_struct(c).learning_repeats))
    xline(mean(corr_vals(c,:)),'LineWidth',1,'Color',colors(c),'DisplayName',sprintf('Mean Correlation After %i Repeats',curr_struct(c).learning_repeats))
end
legend()
xlabel('Correlation Value')
ylabel('Number of Current Sequences')
title(sprintf('Correlation of Sequences to %i Repeat(s) Template \n Before and After STDP',spikes_struct(compare_to).repeat_number))

clear colors c
%% Visualize true vs replay with gaussian conv. of replays
%Max correlation indices and values
[~, max_ind] = max(corr_vals,[],2);
%Convolution parameters
bin_size = 20; %width of Gaussian kernel in bins
x = round(-bin_size/2):round(bin_size/2); %x-values for Gaussian convolution kernel
norm_conv = normpdf(x,0,bin_size/6); %Gaussian convolution kernel
%Figure
figure;
subplot(2,2,1)
spikes_x = unique(spikes_struct(compare_to).spikes_x,'stable');
imagesc(spikes_struct(compare_to).spikes_Vm(spikes_x,1:spikes_struct(1).max_dt))
title('Template Sequence')
ax1 = subplot(2,2,2);
times_1 = curr_struct(1).events(max_ind(1)).times;
seq_spikes_1 = curr_struct(1).spikes_Vm(spikes_x,times_1(1):times_1(2));
seq_1_conv = [];
for i = 1:length(spikes_x)
    bin_conv = conv(norm_conv,seq_spikes_1(i,:));
    seq_1_conv = [seq_1_conv; bin_conv];
end    
clear times_1 seq_spikes_1 i bin_conv
imagesc(seq_1_conv)
title(sprintf('Sequence from %i Repeat',curr_struct(1).learning_repeats))
ax2 = subplot(2,2,3);
times_2 = curr_struct(2).events(max_ind(2)).times;
seq_spikes_2 = curr_struct(2).spikes_Vm(spikes_x,times_2(1):times_2(2));
seq_2_conv = [];
for i = 1:length(spikes_x)
    bin_conv = conv(norm_conv,seq_spikes_1(i,:));
    seq_2_conv = [seq_2_conv; bin_conv];
end    
clear times_2 seq_spikes_2 i bin_conv
imagesc(seq_2_conv)
title(sprintf('Sequence from %i Repeat',curr_struct(2).learning_repeats))
ax3 = subplot(2,2,4);
times_3 = curr_struct(3).events(max_ind(3)).times;
seq_spikes_3 = curr_struct(3).spikes_Vm(spikes_x,times_3(1):times_3(2));
seq_3_conv = [];
for i = 1:length(spikes_x)
    bin_conv = conv(norm_conv,seq_spikes_1(i,:));
    seq_3_conv = [seq_3_conv; bin_conv];
end
clear times_3 seq_spikes_3 i bin_conv
imagesc(seq_3_conv)
title(sprintf('Sequence from %i Repeat',curr_struct(3).learning_repeats))
linkaxes([ax1,ax2,ax3])
clear seq_1_conv seq_2_conv seq_3_conv ax1 ax2 ax3 x bin_size norm_conv max_ind