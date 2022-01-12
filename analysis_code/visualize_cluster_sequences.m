%% Visualize cluster sequences
%This code visualizes the sequence of clusters a particular spike
%sequence progresses through.

%Select and load specific network data to analyze
% net_save_path = uigetdir('/Users/hannahgermaine/Documents/PhD/','Select Network Save Folder'); %Have user input where they'd like the output stored
% slashes = find(net_save_path == '/');
% save_path = net_save_path(1:slashes(end));
% load(strcat(save_path,'/parameters.mat'))
% load(strcat(net_save_path,'/network_cluster_sequences.mat'))
% load(strcat(net_save_path,'/network_var.mat'))

V_m = network_var(1).V_m;
spikes_V_m = V_m >= -50*10^(-3);

%Grab relevant information
[~,inits] = size(network_cluster_sequences); %Grab data sizes

%Visualize cluster sequences
seq_type = {'clusters','movsum','normalized_clusters','normalized_cluster_mov_sum'};

%Create image save path
cluster_save_path = strcat(net_save_path,'/cluster_plots/');
if ~isfolder(cluster_save_path)
    mkdir(cluster_save_path);
end

%Create and Save All Cluster Plots
for s = seq_type
    s_type = s{1};
    for k = 1:inits
        if ~isempty(network_cluster_sequences(k).events)
            [num_events,~] = size(network_cluster_sequences(k).events);
            try %Sometimes event information is saved, but the event falls outside a desireable range so cluster information is not saved
                sequences = network_cluster_sequences(k).clusters(1);
                sequence_names = fieldnames(sequences);
                f = figure;
                axes = [];
                for e_i = 1:num_events
                    event_times = network_cluster_sequences(k).events(e_i,:);
                    event_length = event_times(2) - event_times(1);
                    ax = subplot(1,num_events,e_i);
                    axes(end+1) = ax; %#ok<SAGROW>
                    if strcmp(s_type,'clusters')
                        imagesc(network_cluster_sequences(k).clusters.(sequence_names{e_i}))
                    elseif strcmp(s_type,'movsum')
                        imagesc(network_cluster_sequences(k).movsum.(sequence_names{e_i}))
                    elseif strcmp(s_type,'normalized_clusters')
                        imagesc(network_cluster_sequences(k).normalized_clusters.(sequence_names{e_i}))
                    elseif strcmp(s_type,'normalized_cluster_mov_sum')
                        imagesc(network_cluster_sequences(k).normalized_cluster_mov_sum.(sequence_names{e_i}))
                    end
                    xticks(round(linspace(1,event_length,20))) %20 ticks will be displayed
                    xt = get(gca,'XTick');
                    xtlbl = round(linspace(event_times(1)*parameters.dt,event_times(2)*parameters.dt,numel(xt)),2);
                    xlabel('Time (s)')
                    ylabel('Initialization Number')    
                    title(strcat('Event #',string(e_i)))
                    set(gca, 'XTick',xt, 'XTickLabel',xtlbl)
                end
                if strcmp(s_type,'clusters')
                    sgtitle(strcat('Cluster Sequence for Initialization #',string(k)))
                elseif strcmp(s_type,'movsum')
                    sgtitle(strcat('Moving Sum Cluster Sequence for Initialization #',string(k)))
                elseif strcmp(s_type,'normalized_clusters')
                    sgtitle(strcat('Normalized Cluster Sequence for Initialization #',string(k)))
                elseif strcmp(s_type,'normalized_cluster_mov_sum')
                    sgtitle(strcat('Normalized Moving Sum Cluster Sequence for Initialization #',string(k)))
                end
                %Save Figure
                savefig(f,strcat(cluster_save_path,'/init_',string(k),'_',string(s_type),'_cluster_sequence.fig'))
                saveas(f,strcat(cluster_save_path,'/init_',string(k),'_',string(s_type),'_cluster_sequence.jpg'))
                close(f)
            end
        end
    end
end

%% Visualize Re-ordered Cluster Sequences
%If we look at the sequence of clusters re-ordered by first appearance and
%by maximal representation in the movsum, do we visually see a pattern?

