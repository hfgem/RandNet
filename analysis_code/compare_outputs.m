%Code to compare results from neuron and current initializations

%% Set up variables

%Get paths of data
msgbox('Select folder where neuron initialization results are stored.')
neuron_init_path = uigetdir('/Users/hannahgermaine/Documents/PhD/');
msgbox('Select folder where current initialization results are stored.')
current_init_path = uigetdir('/Users/hannahgermaine/Documents/PhD/');

%Get number and name of network folders per init type
neuron_network_list = {dir(strcat(neuron_init_path,'/network*')).name};
current_network_list = {dir(strcat(current_init_path,'/network*')).name};

%% Pull Data - Spike Sequences
%Here we pull all of the sequences generated into easy-to-work-with
%variables

%Create storage variables
neuron_event_lengths = [];
current_event_lengths = [];
neuron_spike_sequences = struct;
current_spike_sequences = struct;

%First pull event length data on neuron initialization
track_i = 1;
for i = 1:length(neuron_network_list)
    data_folder = strcat(neuron_init_path,'/',neuron_network_list{i});
    try
        load(strcat(data_folder,'/network_spike_sequences.mat'))
        for j= 1:length(network_spike_sequences)
            if ~isempty(network_spike_sequences(j).spike_order)
                try
                    %First store event lengths
                    neuron_event_lengths = [neuron_event_lengths, network_spike_sequences(j).event_lengths]; %#ok<SAGROW>
                catch
                    disp('No Event Lengths')
                end
                try
                    %Next store spike sequences
                    neuron_spike_sequences(track_i).spike_order = network_spike_sequences(j).spike_order;
                    %Next store spike ranks
                    neuron_spike_sequences(track_i).spike_ranks = network_spike_sequences(j).spike_ranks;
                    %Finally store nonspiking neurons
                    neuron_spike_sequences(track_i).nonspiking_neurons = network_spike_sequences(j).nonspiking_neurons;
                catch
                    disp('No spike order/rank/nonspiking_neurons')
                end
                track_i = track_i + 1;
            end    
        end
        clear network_spike_sequences
    catch
        disp("No Spike Sequence File")
    end
end  

track_i = 1;
%Second pull event length data on current initialization
for i = 1:length(neuron_network_list)
    data_folder = strcat(current_init_path,'/',neuron_network_list{i});
    try
        load(strcat(data_folder,'/network_spike_sequences.mat'))
        for j= 1:length(network_spike_sequences)
            if ~isempty(network_spike_sequences(j).spike_order)
                try
                    %First store event lengths
                    current_event_lengths = [current_event_lengths, network_spike_sequences(j).event_lengths]; %#ok<SAGROW>
                catch
                    disp('No Event Lengths')
                end
                try
                    %Next store spike sequences
                    current_spike_sequences(track_i).spike_order = network_spike_sequences(j).spike_order;
                    %Next store spike ranks
                    current_spike_sequences(track_i).spike_ranks = network_spike_sequences(j).spike_ranks;
                    %Finally store nonspiking neurons
                    current_spike_sequences(track_i).nonspiking_neurons = network_spike_sequences(j).nonspiking_neurons;
                catch
                    disp('No spike order/rank/nonspiking_neurons')
                end
                track_i = track_i + 1;
            end    
        end
        clear network_spike_sequences
    catch
        disp("No Spike Sequence File")
    end
end

clear track_i i data_folder j network_spike_sequences

save(strcat(neuron_init_path,'/neuron_event_lengths.mat'),'neuron_event_lengths')
save(strcat(current_init_path,'/current_event_lengths.mat'),'current_event_lengths')
save(strcat(neuron_init_path,'/neuron_spike_sequences.mat'),'neuron_spike_sequences')
save(strcat(current_init_path,'/current_spike_sequences.mat'),'current_spike_sequences')

%% Visualize Event Lengths

%Load data if not already in workspace
load(strcat(neuron_init_path,'/neuron_event_lengths.mat'))
load(strcat(current_init_path,'/current_event_lengths.mat'))

%Select save folder for comparison results
msgbox('Select folder where results are stored.')
save_path = uigetdir('/Users/hannahgermaine/Documents/PhD/');

%First compare the event lengths
f = figure;
h1 = histogram(neuron_event_lengths,10);
xline(mean(neuron_event_lengths),'r--')
hold on
h2 = histogram(current_event_lengths,20);
xline(mean(current_event_lengths),'b--')
legend('Neuron Initialization','Mean Neuron Init. Event Length','Current Initialization','Mean Current Init. Event Length')
title('Histograms of Event Lengths')
xlabel('Event Length in Seconds')
ylabel('Number of Events')
savefig(f,strcat(save_path,'/','comparing_neuron_current_best_nets.fig'))
saveas(f,strcat(save_path,'/','comparing_neuron_current_best_nets.jpg'))
saveas(f,strcat(save_path,'/','comparing_neuron_current_best_nets.svg'))

%% Generate Shuffled Spike Sequences for Comparison

%Load datasets if not already in workspace
load(strcat(neuron_init_path,'/neuron_spike_sequences.mat'))
load(strcat(current_init_path,'/current_spike_sequences.mat'))

%Load parameters for network size info
load(strcat(neuron_init_path,'/parameters.mat'))
n_neur = parameters.n;
clear parameters
load(strcat(current_init_path,'/parameters.mat'))
n_curr = parameters.n;
clear parameters

%Set parameters
shuffle_n = 100;
neuron_viable_inits = 1:length(neuron_spike_sequences);
current_viable_inits = 1:length(current_spike_sequences);

%Generate shuffled data
[neuron_shuffled_spike_sequences] = generate_shuffled_trajectories(n_neur,...
    shuffle_n, neuron_spike_sequences, neuron_viable_inits);
[current_shuffled_spike_sequences] = generate_shuffled_trajectories(n_curr,...
    shuffle_n, current_spike_sequences, current_viable_inits);
save(strcat(neuron_init_path,'/neuron_shuffled_spike_sequences.mat'),'neuron_shuffled_spike_sequences','-v7.3')
save(strcat(current_init_path,'/current_shuffled_spike_sequences.mat'),'current_shuffled_spike_sequences','-v7.3')

%% Calculate Spike Sequence Correlations

%Calculate trajectory similarity within the neuron results and within the
%current results and store

% %Load datasets if not already in workspace
% load(strcat(neuron_init_path,'/neuron_spike_sequences.mat'))
% load(strcat(current_init_path,'/current_spike_sequences.mat'))
% load(strcat(neuron_init_path,'/neuron_shuffled_spike_sequences.mat'))
% load(strcat(current_init_path,'/current_shuffled_spike_sequences.mat'))
% 
% %Load parameters for network size info
% load(strcat(neuron_init_path,'/parameters.mat'))
% n_neur = parameters.n;
% clear parameters
% load(strcat(current_init_path,'/parameters.mat'))
% n_curr = parameters.n;
% clear parameters
% 
% %Select save folder for comparison results
% msgbox('Select folder where results are stored.')
% save_path = uigetdir('/Users/hannahgermaine/Documents/PhD/');

%Set parameters
shuffle_n = length(neuron_shuffled_spike_sequences);
shuffled_viable_inits = [1:shuffle_n];
neuron_viable_inits = 1:length(neuron_spike_sequences);
current_viable_inits = 1:length(current_spike_sequences);

%   Spearman's Rank Correlations

%____Neuron Data Ranks____
[n_ranks, n_ranks_mod] = calculate_trajectory_similarity_spearmans(n_neur, ...
    neuron_viable_inits, neuron_spike_sequences);
    %Remove NaN Values
n_ranks_mod(isnan(n_ranks_mod)) = 0;
    %Turn Rank Matrices Into Vectors of Unique Ranks
n_ranks_vec = nonzeros(triu(n_ranks,1)); %all above diagonal
n_ranks_mod_vec = nonzeros(triu(n_ranks_mod,1)); %all above diagonal
save(strcat(neuron_init_path,'/n_ranks.mat'),'n_ranks','-v7.3')
save(strcat(neuron_init_path,'/n_ranks_mod.mat'),'n_ranks_mod','-v7.3')


%____Current Data Ranks____
[c_ranks, c_ranks_mod] = calculate_trajectory_similarity_spearmans(n_curr, ...
    current_viable_inits, current_spike_sequences);
    %Remove NaN Values
c_ranks_mod(isnan(c_ranks_mod)) = 0;
    %Turn Rank Matrices Into Vectors of Unique Ranks
c_ranks_vec = nonzeros(triu(c_ranks,1)); %all above diagonal
c_ranks_mod_vec = nonzeros(triu(c_ranks_mod,1)); %all above diagonal
save(strcat(current_init_path,'/c_ranks.mat'),'c_ranks','-v7.3')
save(strcat(current_init_path,'/c_ranks_mod.mat'),'c_ranks_mod','-v7.3')

%____Shuffled Data Ranks____
%   Get neuron shuffled ranks
[n_shuffled_ranks, n_shuffled_ranks_mod] = calculate_trajectory_similarity_spearmans(n_neur, ...
    shuffled_viable_inits, neuron_shuffled_spike_sequences);
n_shuffled_ranks_vec = nonzeros(triu(n_shuffled_ranks,1)); %all above diagonal
n_shuffled_ranks_mod_vec = nonzeros(triu(n_shuffled_ranks_mod,1)); %all above diagonal
save(strcat(neuron_init_path,'/n_shuffled_ranks.mat'),'n_shuffled_ranks','-v7.3')
save(strcat(neuron_init_path,'/n_shuffled_ranks_mod.mat'),'n_shuffled_ranks_mod','-v7.3')

%   Get current shuffled ranks
[c_shuffled_ranks, c_shuffled_ranks_mod] = calculate_trajectory_similarity_spearmans(n_curr, ...
    shuffled_viable_inits, current_shuffled_spike_sequences);
c_shuffled_ranks_vec = nonzeros(triu(c_shuffled_ranks,1)); %all above diagonal
c_shuffled_ranks_mod_vec = nonzeros(triu(c_shuffled_ranks_mod,1)); %all above diagonal
save(strcat(current_init_path,'/c_shuffled_ranks.mat'),'c_shuffled_ranks','-v7.3')
save(strcat(current_init_path,'/c_shuffled_ranks_mod.mat'),'c_shuffled_ranks_mod','-v7.3')

%%  Visualize Sequence Correlations

load(strcat(neuron_init_path,'/n_ranks.mat'))
load(strcat(neuron_init_path,'/n_ranks_mod.mat'))
load(strcat(current_init_path,'/c_ranks.mat'))
load(strcat(current_init_path,'/c_ranks_mod.mat'))
load(strcat(neuron_init_path,'/n_shuffled_ranks.mat'))
load(strcat(neuron_init_path,'/n_shuffled_ranks_mod.mat'))
load(strcat(current_init_path,'/c_shuffled_ranks.mat'))
load(strcat(current_init_path,'/c_shuffled_ranks_mod.mat'))

%Turn Rank Matrices Into Vectors of Unique Ranks
%   Grab only values above the diagonal
n_ranks_vec = nonzeros(triu(n_ranks,1));
n_ranks_mod_vec = nonzeros(triu(n_ranks_mod,1));
c_ranks_vec = nonzeros(triu(c_ranks,1));
c_ranks_mod_vec = nonzeros(triu(c_ranks_mod,1));
n_shuffled_ranks_vec = nonzeros(triu(n_shuffled_ranks,1));
n_shuffled_ranks_mod_vec = nonzeros(triu(n_shuffled_ranks_mod,1));
c_shuffled_ranks_vec = nonzeros(triu(c_shuffled_ranks,1));
c_shuffled_ranks_mod_vec = nonzeros(triu(c_shuffled_ranks_mod,1));

%_____Plot Histograms with Percentiles_____
f = figure;
ax1 = subplot(2,2,1);
histogram(n_shuffled_ranks_vec,'DisplayName','Shuffled Values')
hold on
histogram(n_ranks_vec,'DisplayName','Real Data Values')
legend()
title({'Neuron Initialization SRC Rhos','with Nonspiking Neurons at End'})
ax2 = subplot(2,2,2);
histogram(n_shuffled_ranks_mod_vec,'DisplayName','Shuffled Values')
hold on
histogram(n_ranks_mod_vec,'DisplayName','Real Data Values')
legend()
title({'Neuron Initialization SRC Rhos','Excluding Nonspiking Neurons'})
ax3 = subplot(2,2,3);
histogram(c_shuffled_ranks_vec,'DisplayName','Shuffled Values')
hold on
histogram(c_ranks_vec,'DisplayName','Real Data Values')
legend()
title({'Current Initialization SRC Rhos','with Nonspiking Neurons at End'})
ax4 = subplot(2,2,4);
histogram(c_shuffled_ranks_mod_vec,'DisplayName','Shuffled Values')
hold on
histogram(c_ranks_mod_vec,'DisplayName','Real Data Values')
legend()
title({'Current Initialization SRC Rhos','Excluding Nonspiking Neurons'})
f.Position = [187,387,1112,410];
savefig(f,strcat(save_path,'/','spearmans_rank_percentiles.fig'))
saveas(f,strcat(save_path,'/','spearmans_rank_percentiles.jpg'))
saveas(f,strcat(save_path,'/','spearmans_rank_percentiles.svg'))
close(f)

%% Calculate Statistics Comparing Real to Null

load(strcat(neuron_init_path,'/n_ranks.mat'))
load(strcat(neuron_init_path,'/n_ranks_mod.mat'))
load(strcat(current_init_path,'/c_ranks.mat'))
load(strcat(current_init_path,'/c_ranks_mod.mat'))
load(strcat(neuron_init_path,'/n_shuffled_ranks.mat'))
load(strcat(neuron_init_path,'/n_shuffled_ranks_mod.mat'))
load(strcat(current_init_path,'/c_shuffled_ranks.mat'))
load(strcat(current_init_path,'/c_shuffled_ranks_mod.mat'))

%Turn Rank Matrices Into Vectors of Unique Ranks
%   Grab only values above the diagonal
n_ranks_vec = nonzeros(triu(n_ranks,1));
n_ranks_mod_vec = nonzeros(triu(n_ranks_mod,1));
c_ranks_vec = nonzeros(triu(c_ranks,1));
c_ranks_mod_vec = nonzeros(triu(c_ranks_mod,1));
n_shuffled_ranks_vec = nonzeros(triu(n_shuffled_ranks,1));
n_shuffled_ranks_mod_vec = nonzeros(triu(n_shuffled_ranks_mod,1));
c_shuffled_ranks_vec = nonzeros(triu(c_shuffled_ranks,1));
c_shuffled_ranks_mod_vec = nonzeros(triu(c_shuffled_ranks_mod,1));

%____Calculate Similarity to Null Distributions Via Kolmogorov-Smirnov____
%____2-sample test____
[h_n,p_n] = kstest2(n_ranks_vec,n_shuffled_ranks_vec);
[h_c,p_c] = kstest2(c_ranks_vec,c_shuffled_ranks_vec);
[h_n_mod,p_n_mod] = kstest2(n_ranks_mod_vec,n_shuffled_ranks_mod_vec);
[h_c_mod,p_c_mod] = kstest2(c_ranks_mod_vec,c_shuffled_ranks_mod_vec);

%____Calculate Which True Values are Above a Percentile in the____
%____Shuffled Dataset____
%First calculate percentile cutoffs
n_prct_90 = prctile(n_shuffled_ranks_vec,90);
n_prct_95 = prctile(n_shuffled_ranks_vec,95);
n_prct_90_mod = prctile(n_shuffled_ranks_mod_vec,90);
n_prct_95_mod = prctile(n_shuffled_ranks_mod_vec,95);
c_prct_90 = prctile(c_shuffled_ranks_vec,90);
c_prct_95 = prctile(c_shuffled_ranks_vec,95);
c_prct_90_mod = prctile(c_shuffled_ranks_mod_vec,90);
c_prct_95_mod = prctile(c_shuffled_ranks_mod_vec,95);
%Indices of Results Above Percentile Values
n_ranks_90 = n_ranks_vec > n_prct_90;
n_ranks_95 = n_ranks_vec > n_prct_95;
n_ranks_90_mod = n_ranks_mod_vec > n_prct_90_mod;
n_ranks_95_mod = n_ranks_mod_vec > n_prct_95_mod;
c_ranks_90 = c_ranks_vec > c_prct_90;
c_ranks_95 = c_ranks_vec > c_prct_95;
c_ranks_90_mod = c_ranks_mod_vec > c_prct_90_mod;
c_ranks_95_mod = c_ranks_mod_vec > c_prct_95_mod;
%Fractions of Results Above Percentile Values
n_frac_90 = sum(n_ranks_90)/length(n_ranks_90);
n_frac_95 = sum(n_ranks_95)/length(n_ranks_95);
n_frac_90_mod = sum(n_ranks_90_mod)/length(n_ranks_90_mod);
n_frac_95_mod = sum(n_ranks_95_mod)/length(n_ranks_95_mod);
c_frac_90 = sum(c_ranks_90)/length(c_ranks_90);
c_frac_95 = sum(c_ranks_95)/length(c_ranks_95);
c_frac_90_mod = sum(c_ranks_90_mod)/length(c_ranks_90_mod);
c_frac_95_mod = sum(c_ranks_95_mod)/length(c_ranks_95_mod);