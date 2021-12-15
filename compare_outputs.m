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

%Create storage variables
neuron_event_lengths = [];
current_event_lengths = [];
neuron_spike_sequences = struct;
current_spike_sequences = struct;
%% Pull Data

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
                    neuron_event_lengths(end+1) = network_spike_sequences(j).event_lengths; %#ok<SAGROW>
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
                    current_event_lengths(end+1) = network_spike_sequences(j).event_lengths; %#ok<SAGROW>
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

%% Calculate Sequence Comparisons

%% Visualize Comparisons

%First compare the event lengths
figure;
h1 = histogram(neuron_event_lengths);
xline(mean(neuron_event_lengths),'r--')
hold on
h2 = histogram(current_event_lengths);
xline(mean(current_event_lengths),'b--')
legend('Neuron Initialization','Mean Neuron Init. Event Length','Current Initialization','Mean Current Init. Event Length')
title('Histograms of Event Lengths')
xlabel('Event Length in Seconds')
ylabel('Number of Events')






