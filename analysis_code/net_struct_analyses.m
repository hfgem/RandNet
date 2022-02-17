%Network Structure Analyses

saveFlag = 0; % 1 to save analysis results

%% Select data to visualize
net_save_path = uigetdir('/Users/hannahgermaine/Documents/PhD/','Select Network Save Folder'); %Have user input where they'd like the output stored
load(strcat(net_save_path,'/network.mat'))
slashes = find(net_save_path == '/');
save_path = net_save_path(1:slashes(end));
load(strcat(save_path,'/parameters.mat'))
scriptFolder = '/netStruct'; % sub-folder so save analysis results to

if saveFlag & ~isfolder([net_save_path, scriptFolder])
    mkdir([net_save_path, scriptFolder]);
end

cluster_mat = network.cluster_mat;
clusters = parameters.clusters;
n = parameters.n;
conns = network.conns;

%% Visualize network structure

%First create indices for each neuron in each cluster - staggering
%locations
for i = 1:clusters
    for j = 1:n
        if i <= clusters/2
            x = i;
            y_base = (clusters/2) - abs((clusters/4) - i);
        else
            x = clusters - i;
            y_base = abs((3*clusters/4) - i);
        end
        cluster_plot_indices(i,j,1) = x + 0.1*randn();
        cluster_plot_indices(i,j,2) = y_base + 0.1*randn();
    end
end

%Next find overlapping neuron scatter positions to draw connecting lines
cluster_plot_connections = [];
for i = 1:n
    clusters_in = find(cluster_mat(:,i));
    if length(clusters_in) > 2
        possible_pairs = nchoosek(clusters_in,2);
        for j = 1:length(possible_pairs)
            index_pair = zeros(2,2);
            index_pair(1,:) = cluster_plot_indices(possible_pairs(j,1),i,:);
            index_pair(2,:) = cluster_plot_indices(possible_pairs(j,2),i,:);
            cluster_plot_connections = cat(3,cluster_plot_connections,index_pair); %#ok<AGROW>
        end    
    elseif length(clusters_in) == 2
        index_pair = zeros(2,2);
        index_pair(1,:) = cluster_plot_indices(clusters_in(1),i,:);
        index_pair(2,:) = cluster_plot_indices(clusters_in(2),i,:);
    end
end
[~, ~, cluster_connections_l] = size(cluster_plot_connections);

%Set colors for each cluster
color_options = jet(clusters);

%Plot
f = figure;
leg_items = [];
leg_names = [];
hold on
for i = 1:clusters
    scat = scatter(cluster_plot_indices(i,:,1),cluster_plot_indices(i,:,2),[],color_options(i,:),'filled');
    leg_items(end+1) = scat;
    leg_names = [leg_names, strcat('Cluster #',string(i))]; %#ok<AGROW>
end
for i = 1:cluster_connections_l
    line(cluster_plot_connections(:,1,i),cluster_plot_connections(:,2,i))
end
legend(leg_items,leg_names)
title('Network Diagram','FontSize',16)
set(gca,'xtick',[],'ytick',[])
if saveFlag
    savefig(f,strcat(net_save_path,scriptFolder,'/','network.fig'))
    saveas(f,strcat(net_save_path,scriptFolder,'/','network.jpg'))
    saveas(f,strcat(net_save_path,scriptFolder,'/','network.svg'))
end

%Plot individual cluster connections
f2 = figure;
hold on
for i = 1:clusters
    within_cluster_connections = [];
    %First find and store intra-cluster connections
    neur_ind = find(cluster_mat(i,:));
    for j = neur_ind %row
        for k = neur_ind %column
            if conns(j,k) == 1
                index_pair = zeros(2,2);
                index_pair(1,:) = cluster_plot_indices(i,j,:);
                index_pair(2,:) = cluster_plot_indices(i,k,:);
                within_cluster_connections = cat(3,within_cluster_connections,index_pair); %#ok<AGROW>
            end    
        end
    end
    [~, ~, cluster_connections_l] = size(within_cluster_connections);
    %Plot
    scat = scatter(cluster_plot_indices(i,:,1),cluster_plot_indices(i,:,2),[],color_options(i,:),'filled');
    for i = 1:cluster_connections_l
        line(within_cluster_connections(:,1,i),within_cluster_connections(:,2,i))
    end
end 
title('Intra-Cluster Connectivity','FontSize',16)
set(gca,'xtick',[],'ytick',[])
if saveFlag
    savefig(f2,strcat(net_save_path,scriptFolder,'/','network_clusters.fig'))
end

%% Creat Plots of Network Properties
%{
%Select data to visualize
net_save_path = uigetdir('/Users/hannahgermaine/Documents/PhD/','Select Network Save Folder'); %Have user input where they'd like the output stored
load(strcat(net_save_path,'/network.mat'))
slashes = find(net_save_path == '/');
save_path = net_save_path(1:slashes(end));
load(strcat(save_path,'/parameters.mat'))

cluster_mat = network.cluster_mat;
clusters = parameters.clusters;
n = parameters.n;
conns = network.conns;
%}

%Plot Cluster Participation as Imagesc
f = figure;
imagesc(cluster_mat)
colorbar('Ticks',[0,1],'TickLabels',{'No','Yes'})
colormap(summer)
xlabel('Neuron Index')
ylabel('Cluster Index')
title('Cluster Participation')
f.Position = [440,720,866,85];
if saveFlag
    savefig(f,strcat(net_save_path,scriptFolder,'/','cluster_participation_visual.fig'))
    saveas(f,strcat(net_save_path,scriptFolder,'/','cluster_participation_visual.jpg'))
    saveas(f,strcat(net_save_path,scriptFolder,'/','cluster_participation_visual.svg'))
    %close(f)   
end

%Plot Histogram of number of clusters each neuron participates in
num_clusters = sum(cluster_mat,1);
mean_clusters = mean(num_clusters);
f2 = figure;
histogram(num_clusters,'FaceColor',[1.0000,1.0000,0.4000])
xline(mean_clusters,'label',strcat('Mean = ',string(mean_clusters)),...
    'Color',[0,0.5000,0.4000],'LineWidth',1)
ylabel('Number of Neurons')
xlabel('Number of Clusters')
title('Most Neurons Participate in Multiple Clusters')
if saveFlag
    savefig(f2,strcat(net_save_path,scriptFolder,'/','cluster_participation_histogram.fig'))
    saveas(f2,strcat(net_save_path,scriptFolder,'/','cluster_participation_histogram.jpg'))
    saveas(f2,strcat(net_save_path,scriptFolder,'/','cluster_participation_histogram.svg'))
    %close(f2)
end

%Plot Histogram of number of neurons that participate in a cluster overlap
pairs = nchoosek(1:clusters,2);
overlap = zeros(1,length(pairs));
for i = 1:length(pairs)
    overlap(1,i) = sum(sum(cluster_mat(pairs(i,:),:),1) == 2); %neurons in both clusters
end
mean_overlap = mean(overlap);
f3 = figure;
histogram(overlap,'FaceColor',[1.0000,0.7812,0.4975])
xline(mean_overlap,'label',strcat('Mean = ',string(mean_overlap)),...
    'Color',[0,0.5000,0.4000],'LineWidth',1)
ylabel('Number of Cluster Pairs')
xlabel('Number of Neurons')
title('Number of Neurons in Cluster Overlap')
if saveFlag
    savefig(f3,strcat(net_save_path,scriptFolder,'/','cluster_overlap_histogram.fig'))
    saveas(f3,strcat(net_save_path,scriptFolder,'/','cluster_overlap_histogram.jpg'))
    saveas(f3,strcat(net_save_path,scriptFolder,'/','cluster_overlap_histogram.svg'))
    %close(f3)
end

clear f f2
%% Creat plot of overlap in 2 clusters

%Select data to visualize
%{
net_save_path = uigetdir('/Users/hannahgermaine/Documents/PhD/','Select Network Save Folder'); %Have user input where they'd like the output stored
load(strcat(net_save_path,'/network.mat'))
slashes = find(net_save_path == '/');
save_path = net_save_path(1:slashes(end));
load(strcat(save_path,'/parameters.mat'))

cluster_mat = network.cluster_mat;
clusters = parameters.clusters;
cluster_n = parameters.cluster_n;
n = parameters.n;
conns = network.conns;
%}

clust_ind = randperm(clusters,2); %Randomly pick 2 clusters for visualization
neur_ind = find(sum(cluster_mat(clust_ind,:),1) == 2); %neurons in both clusters
colors = [[0,0.5000,0.7500];[1.0000,0.7812,0.4975];[0,0.5000,0.4000]];

%Store indices and colors of scatter points
plot_ind = zeros(n,2); %index of each neuron scatter plot
scatter_col = zeros(n,3); %color of each neuron scatter plot

%Assign neuron color and location based on its cluster participation
for i = 1:n
    if sum(neur_ind == i) %is in both clusters
        x_center = 1.5;
        y_center = 1;
        scatter_col(i,:) = colors(3,:);
        plot_ind(i,1) = x_center + 0.2*(rand() - 0.5);
        plot_ind(i,2) = y_center + 0.75*(rand() - 0.5);
    else
        in_clust = find(cluster_mat(clust_ind,i));
        if ~isempty(in_clust) %neuron is in one of the clusters
            if in_clust == 2 %set centers of point clouds
                x_center = 2;
                y_center = 1;
                scatter_col(i,:) = colors(2,:);
            else
                x_center = 1;
                y_center = 1;
                scatter_col(i,:) = colors(1,:);
            end
            plot_ind(i,1) = x_center + 0.5*(rand() - 0.5);
            plot_ind(i,2) = y_center + 0.5*(rand() - 0.5);
        end
    end
end

%Plot Cluster Connectivities
f = figure;
hold on
for i = 1:2
    c_i = clust_ind(i);
    line_color = colors(i,:);
    for j = 1:n
        if cluster_mat(c_i,j) == 1
            for k = 1:n
                if cluster_mat(c_i,k) == 1
                    if logical(conns(j,k)) == 1
                        line([plot_ind(j,1),plot_ind(k,1)],[plot_ind(j,2),...
                            plot_ind(k,2)],'LineStyle',':','Color',line_color)
                    end
                end
            end
        end
    end
end
scatter(plot_ind(:,1),plot_ind(:,2),50,scatter_col,'filled')
set(gca,'xtick',[],'ytick',[])
if saveFlag
    savefig(f,strcat(net_save_path,scriptFolder,'/','cluster_overlap_vis.fig'))
    saveas(f,strcat(net_save_path,scriptFolder,'/','cluster_overlap_vis.jpg'))
    saveas(f,strcat(net_save_path,scriptFolder,'/','cluster_overlap_vis.svg'))
    close(f)
end