%Testing how STDP applied to one sequence affects other sequences

%______ABOUT______:
%This code tests how repeating one sequence through re-initialization
%affects other sequences. The code also makes use of previous network
%structures and parameters that generated successful sequences (from
%lif_network_postrotation.m) to load a constant network structure for
%testing.
%__________________

%% Initialization

%Select folder of good network structure and parameters to use in tests
load_path = uigetdir('/Users/hannahgermaine/Documents/PhD/','Select network.m Load Folder'); %Have user input where they'd like network structure to be loaded from

% %_________________________________
% %___Load Independent Parameters___
% %_________________________________
load(strcat(load_path,'/network.mat'))
slashes = find(load_path == '/');
param_path = load_path(1:slashes(end));
load(strcat(param_path,'/parameters.mat'))

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

%How many tests of different initializations to run
if strcmp(parameters.type,'cluster')
    test_val_max = parameters.clusters; %Ensures every cluster is initialized
else
    test_val_max = 10; %This value can be modified as you see fit
end
%save for easy calculations
parameters.('test_val_max') = test_val_max;

%Adding an input current to all cells (one of the options must be uncommented)
x_in = [0:parameters.dt:parameters.t_max];
% %Rhythmic current input: (uncomment if desired)
% I_coeff = 0; %5.1*10^(-10); %set to 0 for no input current
% I_Hz = 1; %frequency of input - 1 Hz or 60 bpm, 4-8 Hz for theta rhythm
% I_in = I_coeff*(0.5 + 0.5*sin(I_Hz*2*pi*x_in)); %Approximately theta wave input current
%Noisy input: (uncomment if desired)
G_in = parameters.G_coeff*randn(parameters.n,parameters.t_steps+1)*parameters.G_scale;

%save for easy calculations
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

%____________________________________
%___Assign STDP Rules Desired________
%____________________________________
%Here you can modify the STDP rules previously saved
parameters.type = 'neuron'; %Each sequence initialization will be through setting neurons to threshold
parameters.tau_stdp = 1*10^(-3); %STDP time constant (s)
parameters.connectivity_gain = 0.005; %amount to increase or decrease connectivity by with each STDP rule (more at the range of 0.002-0.005)
%Here you set the simulation rules
num_repeat_tests = 10; %how many repeat values to test
max_repeats = 50; %maximum number of repeats to test
num_repeats = 0:max_repeats/num_repeat_tests:max_repeats+1;  %vector of repeats to test
num_repeats(1) = []; %update first value to be 1 repeat rather than 0

%% Test STDP's Effect on Other Sequences

%Run learning on one of the initializations the number of times
%specified in num_repeats and save the final sequence of the
%initialization used, and then run the initializations for other sequences
%and save them

%Select simulation save folder
save_path = uigetdir('/Users/hannahgermaine/Documents/PhD/','Select Save Folder'); %Have user input where they'd like network structure to be loaded from

%Create a structure file to store all the sequence results
STDP_sequences = struct;
STDP_V_m = struct;

init_seed = [2,4,1,8,9,6,10,3,5,7];

for i = 1:num_repeat_tests
    %Create a network copy to modify
    network_copy = network;
    
    %Set up storage
    STDP_sequences(i).num_repeats = num_repeats(i);
    STDP_sequences(i).init_seed = init_seed;
    STDP_sequences(i).spike_order = zeros(parameters.n,parameters.test_val_max);
    STDP_sequences(i).nonspiking = zeros(parameters.n,parameters.test_val_max);
    
    %Run for selected number of repeats
    %init_seed(1) will be the one that is repeated
    for j = 1:num_repeats(i)
        %Create Storage Variables
        I_syn = zeros(parameters.n,parameters.t_steps+1); %synaptic current emitted by each neuron at each timestep (A)
        %synaptic conductances for each neuron at each timestep
        G_syn_I = zeros(parameters.n,parameters.t_steps+1); %conductance for presynaptic inhibitory (S)
        G_syn_E = zeros(parameters.n,parameters.t_steps+1); %conductance for presynaptic excitatory (S)
        V_m = zeros(parameters.n,parameters.t_steps+1); %membrane potential for each neuron at each timestep
        V_m(:,1) = parameters.V_reset + randn([parameters.n,1])*(10^(-3)); %set all neurons to baseline reset membrane potential with added noise
        G_sra = zeros(parameters.n,parameters.t_steps+1); %refractory conductance for each neuron at each timestep
        
        %Run simulation
        [V_m, ~, ~, ~, ~, conns_update] = lif_sra_calculator_postrotation(...
            parameters, init_seed(1), network_copy, I_syn, G_syn_I, G_syn_E, V_m, G_sra);
        %update network file with new connectivity matrix from STDP
        network_copy.conns = conns_update;
    end
    STDP_sequences(i).conns = network_copy.conns; %Save the final connectivity matrix
    STDP_V_m(i).(strcat('V_m_',string(1))) = V_m;
    
    %Store sequence from init_seed(1)
    spikes_V_m = V_m >= parameters.V_th;
    [spikes_x,~] = find(spikes_V_m);
    spike_order = unique(spikes_x,'stable');
    nonspiking = sum(spikes_V_m,2) == 0;
    STDP_sequences(i).spike_order(:,1) = [spike_order; find(nonspiking)];
    STDP_sequences(i).nonspiking(:,1) = nonspiking;
    
    %Clear unnecessary variables
    clear j I_syn G_syn_I G_syn_E V_m G_sra conns_update spikes_V_m ...
        spikes_x spike_order nonspiking
    
    %Run other sequence initializations
    for j = 2:parameters.test_val_max
        %Create Storage Variables
        I_syn = zeros(parameters.n,parameters.t_steps+1); %synaptic current emitted by each neuron at each timestep (A)
        %synaptic conductances for each neuron at each timestep
        G_syn_I = zeros(parameters.n,parameters.t_steps+1); %conductance for presynaptic inhibitory (S)
        G_syn_E = zeros(parameters.n,parameters.t_steps+1); %conductance for presynaptic excitatory (S)
        V_m = zeros(parameters.n,parameters.t_steps+1); %membrane potential for each neuron at each timestep
        V_m(:,1) = parameters.V_reset + randn([parameters.n,1])*(10^(-3)); %set all neurons to baseline reset membrane potential with added noise
        G_sra = zeros(parameters.n,parameters.t_steps+1); %refractory conductance for each neuron at each timestep
        
        %Run simulation
        [V_m, ~, ~, ~, ~, ~] = lif_sra_calculator_postrotation(...
            parameters, init_seed(j), network_copy, I_syn, G_syn_I, G_syn_E, V_m, G_sra);
        STDP_V_m(i).(strcat('V_m_',string(j))) = V_m;
        
        %Store sequence from this initialization
        spikes_V_m = V_m >= parameters.V_th;
        [spikes_x,~] = find(spikes_V_m);
        spike_order = unique(spikes_x,'stable');
        nonspiking = sum(spikes_V_m,2) == 0;
        STDP_sequences(i).spike_order(:,j) = [spike_order; find(nonspiking)];
        STDP_sequences(i).nonspiking(:,j) = nonspiking;
        
        %Clear unnecessary variables
        clear j I_syn G_syn_I G_syn_E V_m G_sra conns_update spikes_V_m ...
            spikes_x spike_order nonspiking
    end 
end 

%SAVE DATA
save(strcat(save_path,'/STDP_sequences.mat'),'STDP_sequences','-v7.3')
save(strcat(save_path,'/STDP_V_m.mat'),'STDP_V_m','-v7.3')
save(strcat(save_path,'/network.mat'),'network_copy','-v7.3') %contains the updated connectivity matrix

disp('STDP Tests Done')

%% Analyze the Results of STDP on One Sequence
%Load data
%Select folder of STDP datasets
% load_path = uigetdir('/Users/hannahgermaine/Documents/PhD/','Select Save Folder');
% load(strcat(load_path,'/STDP_sequences.mat'))
% load(strcat(load_path,'/STDP_V_m.mat'))
% slashes = find(load_path == '/');
% param_path = load_path(1:slashes(end));
% load(strcat(param_path,'/parameters.mat'))

% load_path = save_path;

%_____________
% Visualize how the sequence changes over the course of STDP learning
%_____________
%Pull /data
num_repeats = [STDP_sequences.num_repeats];
num_to_visualize = length(num_repeats);
same_sequences = zeros(parameters.n,num_to_visualize);
same_non_spiking = zeros(parameters.n,num_to_visualize);
for i = 1:num_to_visualize
    same_sequences(:,i) = STDP_sequences(i).spike_order(:,1);
    same_non_spiking(:,i) = STDP_sequences(i).nonspiking(:,1);
end
%Set up plot
subplot_x = floor(sqrt(num_to_visualize));
subplot_y = ceil(sqrt(num_to_visualize));
f = figure;
axes = [];
%To order all sequences by a particular firing sequence, uncomment the
%following and select which sequence to order by by changing the y-index
% STDP_order = same_sequences(:,1);
for i = 1:num_to_visualize
    ax = subplot(subplot_x,subplot_y,i);
    V_m_spikes = STDP_V_m(i).('V_m_1') >= parameters.V_th;
    STDP_order = same_sequences(:,i); %Comment out if ordering by line 217
    V_m_spikes_reordered = V_m_spikes(STDP_order,:);
    spike_times = find(sum(V_m_spikes_reordered,1));
    large_spike_intervals = find((spike_times(2:end) - spike_times(1:end-1))>500);
    if isempty(large_spike_intervals)
        last_spike_time = spike_times(end);
    else
        last_spike_time = spike_times(large_spike_intervals(1))+100;
    end
    last_spike_time_s = last_spike_time * parameters.dt * 10^3; %convert to ms
    imagesc(V_m_spikes_reordered(:,1:last_spike_time))
    colormap(flip(gray))
    xticks(round(linspace(1,last_spike_time,5))) %20 ticks will be displayed
    xt = get(gca,'XTick');
    xtlbl = round(linspace(0,last_spike_time_s,numel(xt)),2);
    set(gca, 'XTick',xt, 'XTickLabel',xtlbl)
%     xlabel('Time (ms)')%,'FontSize',16)
%     ylabel('Reordered Neuron Number')%,'FontSize',16)
    title(strcat('Repeats = ',string(num_repeats(i))))
    axes = [axes, ax]; %#ok<AGROW>
end
linkaxes(axes)
set(gcf, 'Position', get(0, 'Screensize'));
sgtitle('Sequence as a Function of Learning')
savefig(f,strcat(load_path,'/sequences_after_learning.fig'))
saveas(f,strcat(load_path,'/sequences_after_learning.jpg'))

%_____________
% Calculate the sequence correlations for different amounts of learning
% This calculation excludes nonspiking neurons
%_____________
ranks = zeros(num_to_visualize,num_to_visualize);
for i = 1:num_to_visualize
    X_i = same_sequences(:,i);
    X_nonspike = find(same_non_spiking(:,i));
    for j = 1:num_to_visualize
        Y_i = same_sequences(:,j);
        Y_nonspike = find(same_non_spiking(:,j));
        nonspiking = unique([X_nonspike',Y_nonspike']);
        %Only look at neurons spiking in both
        X_i = setdiff(X_i, nonspiking', 'stable');
        Y_i = setdiff(Y_i, nonspiking', 'stable');
        new_n = (length(X_i)+length(Y_i))/2;
        if new_n ~= 0
            d = zeros(1,parameters.n);
            for k = 1:parameters.n
                %find the ranks of each neuron
                ix_i = find(X_i == k);
                iy_i = find(Y_i == k);
                if ~isempty(ix_i) && ~isempty(iy_i) %if both are spiking
                    d(1,k) = ix_i - iy_i;
                end
            end
            ranks(i,j) = 1 - (6*sum(d.^2))/(new_n*(new_n^2 - 1));
            ranks(j,i) = ranks(i,j);
            if i == j
                ranks(i,j) = 1;
            end
        end
    end
end    
clear i X_i X_nonspike j Y_i Y_nonspike nonspiking new_n d k ix_i iy_i

ranks_vec = nonzeros(triu(ranks,1));

f1 = figure;
subplot(1,2,1)
histogram(ranks_vec,num_to_visualize)
title('Histogram of Unique Rank Correlation Values')
subplot(1,2,2)
imagesc(ranks)
colorbar()
xticks([1:num_to_visualize])
yticks([1:num_to_visualize])
xt = get(gca,'XTick');
yt = get(gca,'YTick');
set(gca, 'XTick',xt, 'XTickLabel',num_repeats)
set(gca, 'YTick',yt, 'YTickLabel',num_repeats)
set(gcf, 'Position', get(0, 'Screensize'));
title('Visualization of Rank Correlation Pairs')
savefig(f1,strcat(load_path,'/sequences_after_learning_correlation.fig'))
saveas(f1,strcat(load_path,'/sequences_after_learning_correlation.jpg'))

% Visualize histograms of correlation values for different learning amounts

%_____________
% Calculate the sequence correlations for different amounts of learning
% ranks_win looks at rank correlations of different sequences compared to
%   each other when one of them has been strengthened through STDP
% ranks_across looks at rank correlations of the same sequence across
%   different learning of a different sequence
%_____________

%First store sequences across inits
[n,num_inits] = size(STDP_sequences(1).spike_order);
sequences_across = zeros(num_inits,n,num_to_visualize);
sequences_across_nonspiking = zeros(num_inits,n,num_to_visualize);
sequences_across_length = zeros(num_inits,num_to_visualize);
for i = 1:num_inits
    for j = 1:num_to_visualize
        sequences_across(i,:,j) = STDP_sequences(j).spike_order(:,i);
        sequences_across_nonspiking(i,:,j) = STDP_sequences(j).nonspiking(:,i);
        sequences_across_length(i,j) = n - sum(STDP_sequences(j).nonspiking(:,i));
    end
end

%Calculate ranks for sequences within a learning amount
ranks_win = zeros(num_to_visualize,num_inits,num_inits);
for i = 1:num_to_visualize
    diff_sequences = STDP_sequences(i).spike_order;
    diff_sequences_nonspiking = STDP_sequences(i).nonspiking;
    for j = 1:num_inits
        X_i = diff_sequences(:,j);
        X_nonspike = find(diff_sequences_nonspiking(:,j));
        for k = 1:num_inits
            if j == k
                ranks_win(i,j,k) = 1;
            else
                Y_i = diff_sequences(:,k);
                Y_nonspike = find(diff_sequences_nonspiking(:,k));
                nonspiking = unique([X_nonspike',Y_nonspike']);
                %Only look at neurons spiking in both
                X_i = setdiff(X_i, nonspiking', 'stable');
                Y_i = setdiff(Y_i, nonspiking', 'stable');
                new_n = (length(X_i)+length(Y_i))/2;
                d = zeros(1,parameters.n);
                for l = 1:parameters.n
                    %find the ranks of each neuron
                    ix_i = find(X_i == l);
                    iy_i = find(Y_i == l);
                    if ~isempty(ix_i) && ~isempty(iy_i) %if both are spiking
                        d(1,l) = ix_i - iy_i;
                    end
                end
                ranks_win(i,j,k) = 1 - (6*sum(d.^2))/(new_n*(new_n^2 - 1));
                ranks_win(i,k,j) = ranks_win(i,j,k);
            end
        end
    end
end
clear i j diff_sequences diff_sequences_nonspiking X_i X_nonspike k Y_i ...
    Y_nonspike nonspiking new_n d l ix_i iy_i

ranks_win_vec = struct;
for i = 1:num_to_visualize
    vec = [];
    for j = 1:num_inits-1
        for k = j+1:num_inits
            vec(end+1) = ranks_win(i,j,k); %#ok<SAGROW>
        end
    end
    ranks_win_vec(i).ranks = vec(~isnan(vec));
    ranks_win_vec(i).avg = mean(ranks_win_vec(i).ranks);
end
clear i vec j k

% Visualize histograms of w/in trial correlations
f2 = figure;
axes = [];
%To order all sequences by a particular firing sequence, uncomment the
%following and select which sequence to order by by changing the y-index
% STDP_order = same_sequences(:,1);
subplot_x = floor(sqrt(num_to_visualize));
subplot_y = ceil(sqrt(num_to_visualize));
for i = 1:num_to_visualize
    ax = subplot(subplot_x,subplot_y,i);
    histogram(ranks_win_vec(i).ranks)
    hold on
    xline(ranks_win_vec(i).avg)
    xlabel('Spearman Rank Correlation')%,'FontSize',16)
    ylabel('Number')%,'FontSize',16)
    title(strcat('Repeats = ',string(num_repeats(i))))
    axes = [axes, ax]; %#ok<AGROW>
end
linkaxes(axes)
set(gcf, 'Position', get(0, 'Screensize'));
sgtitle('Across-Initialization Correlation as a Function of Learning')
savefig(f2,strcat(load_path,'/seq_corr_hist_after_learning.fig'))
saveas(f2,strcat(load_path,'/seq_corr_hist_after_learning.jpg'))
clear i ax axes

% Visualize images of w/in trial correlations
f3 = figure;
axes = [];
%To order all sequences by a particular firing sequence, uncomment the
%following and select which sequence to order by by changing the y-index
subplot_x = floor(sqrt(num_to_visualize));
subplot_y = ceil(sqrt(num_to_visualize));
for i = 1:num_to_visualize
    ax = subplot(subplot_x,subplot_y,i);
    imagesc(squeeze(ranks_win(i,:,:)))
    colorbar()
    xlabel('Number of Initializations')%,'FontSize',16)
    ylabel('Number of Initializations')%,'FontSize',16)
    title(strcat('Repeats = ',string(num_repeats(i))))
    axes = [axes, ax]; %#ok<AGROW>
end
linkaxes(axes)
set(gcf, 'Position', get(0, 'Screensize'));
sgtitle('Across-Initialization Correlation as a Function of Learning')
savefig(f3,strcat(load_path,'/seq_corr_im_after_learning.fig'))
saveas(f3,strcat(load_path,'/seq_corr_im_after_learning.jpg'))
clear i ax axes

%Calculate ranks for sequences across learning rules
ranks_across = zeros(num_inits,num_to_visualize,num_to_visualize);
for i = 1:num_inits
    diff_sequences = squeeze(sequences_across(i,:,:));
    diff_sequences_nonspiking = squeeze(sequences_across_nonspiking(i,:,:));
    for j = 1:num_to_visualize
        X_i = diff_sequences(:,j);
        X_nonspike = find(diff_sequences_nonspiking(:,j));
        for k = 1:num_to_visualize
            if j == k
                ranks_across(i,k,j) = 1;
            else
                Y_i = diff_sequences(:,k);
                Y_nonspike = find(diff_sequences_nonspiking(:,k));
                nonspiking = unique([X_nonspike',Y_nonspike']);
                %Only look at neurons spiking in both
                X_i = setdiff(X_i, nonspiking', 'stable');
                Y_i = setdiff(Y_i, nonspiking', 'stable');
                new_n = (length(X_i)+length(Y_i))/2;
                d = zeros(1,parameters.n);
                for l = 1:parameters.n
                    %find the ranks of each neuron
                    ix_i = find(X_i == l);
                    iy_i = find(Y_i == l);
                    if ~isempty(ix_i) && ~isempty(iy_i) %if both are spiking
                        d(1,l) = ix_i - iy_i;
                    end
                end
                ranks_across(i,j,k) = 1 - (6*sum(d.^2))/(new_n*(new_n^2 - 1));
                ranks_across(i,k,j) = ranks_across(i,j,k) ;
            end
        end
    end
end
clear i j diff_sequences diff_sequences_nonspiking X_i X_nonspike k Y_i ...
    Y_nonspike nonspiking new_n d l ix_i iy_i

ranks_across_vec = struct;
for i = 1:num_to_visualize
    vec = [];
    for j = 1:num_inits-1
        for k = j+1:num_inits
            vec(end+1) = ranks_win(i,j,k); %#ok<SAGROW>
        end
    end
    ranks_across_vec(i).ranks = vec(~isnan(vec));
    ranks_across_vec(i).avg = mean(ranks_across_vec(i).ranks);
end
clear i vec j k

% Visualize images of across trial correlations
f4 = figure;
axes = [];
%To order all sequences by a particular firing sequence, uncomment the
%following and select which sequence to order by by changing the y-index
subplot_x = floor(sqrt(num_inits));
subplot_y = ceil(sqrt(num_inits));
for i = 1:num_inits
    ax = subplot(subplot_x,subplot_y,i);
    imagesc(squeeze(ranks_across(i,:,:)))
    colorbar()
    xticks(1:num_to_visualize)
    yticks(1:num_to_visualize)
    xt = get(gca,'XTick');
    yt = get(gca,'YTick');
    set(gca, 'XTick',xt, 'XTickLabel',num_repeats)
    set(gca, 'YTick',yt, 'YTickLabel',num_repeats)
    xlabel('Number of Repeats')%,'FontSize',16)
    ylabel('Number of Repeats')%,'FontSize',16)
    title([strcat('Initialization = ',string(i));strcat('Sequence length =',string(sequences_across_length(i,1)))])
    axes = [axes, ax]; %#ok<AGROW>
end
linkaxes(axes)
set(gcf, 'Position', get(0, 'Screensize'));
sgtitle('Same Initialization Correlation as a Function of Learning')
savefig(f4,strcat(load_path,'/same_seq_corr_after_learning.fig'))
saveas(f4,strcat(load_path,'/same_seq_corr_after_learning.jpg'))
clear i ax axes

% Sequence length as a function of learning
f5 = figure;
axes = [];
subplot_x = floor(sqrt(num_inits));
subplot_y = ceil(sqrt(num_inits));
for i = 1:num_inits
    ax = subplot(subplot_x,subplot_y,i);
    plot(sequences_across_length(i,:))
    xticks(1:num_to_visualize)
    yticks(1:num_to_visualize)
    xt = get(gca,'XTick');
    yt = get(gca,'YTick');
    set(gca, 'XTick',xt, 'XTickLabel',num_repeats)
    set(gca, 'YTick',yt, 'YTickLabel',num_repeats)
    xlabel('Number of Repeats')
    ylabel('Sequence Length')
    title(strcat('Initialization = ',string(i)))
    axes = [axes, ax]; %#ok<AGROW>
end
set(gcf, 'Position', get(0, 'Screensize'));
sgtitle('Same Initialization Length as a Function of Learning')
savefig(f5,strcat(load_path,'/same_seq_len_after_learning.fig'))
saveas(f5,strcat(load_path,'/same_seq_len_after_learning.jpg'))
clear i ax axes

%%
to_visualize = 2;
str_vis = strcat('V_m_',string(to_visualize));
num_repeats = [STDP_sequences.num_repeats];
num_to_visualize = length(num_repeats);
same_sequences = zeros(parameters.n,num_to_visualize);
same_non_spiking = zeros(parameters.n,num_to_visualize);
for i = 1:num_to_visualize
    same_sequences(:,i) = STDP_sequences(i).spike_order(:,to_visualize);
    same_non_spiking(:,i) = STDP_sequences(i).nonspiking(:,to_visualize);
end
%Set up plot
subplot_x = floor(sqrt(num_to_visualize));
subplot_y = ceil(sqrt(num_to_visualize));
f = figure;
axes = [];
%To order all sequences by a particular firing sequence, uncomment the
%following and select which sequence to order by by changing the y-index
% STDP_order = same_sequences(:,1);
for i = 1:num_to_visualize
    ax = subplot(subplot_x,subplot_y,i);
    V_m_spikes = STDP_V_m(i).(str_vis) >= parameters.V_th;
    STDP_order = same_sequences(:,i); %Comment out if ordering by line 217
    V_m_spikes_reordered = V_m_spikes(STDP_order,:);
    spike_times = find(sum(V_m_spikes_reordered,1));
    large_spike_intervals = find((spike_times(2:end) - spike_times(1:end-1))>500);
    if isempty(large_spike_intervals)
        last_spike_time = spike_times(end);
    else
        last_spike_time = spike_times(large_spike_intervals(1))+100;
    end
    last_spike_time_s = last_spike_time * parameters.dt * 10^3; %convert to ms
    imagesc(V_m_spikes_reordered(:,1:last_spike_time))
    colormap(flip(gray))
    xticks(round(linspace(1,last_spike_time,5))) %20 ticks will be displayed
    xt = get(gca,'XTick');
    xtlbl = round(linspace(0,last_spike_time_s,numel(xt)),2);
    set(gca, 'XTick',xt, 'XTickLabel',xtlbl)
%     xlabel('Time (ms)')%,'FontSize',16)
%     ylabel('Reordered Neuron Number')%,'FontSize',16)
    title(strcat('Repeats = ',string(num_repeats(i))))
    axes = [axes, ax]; %#ok<AGROW>
end
linkaxes(axes)
sgtitle('Sequence as a Function of Learning')