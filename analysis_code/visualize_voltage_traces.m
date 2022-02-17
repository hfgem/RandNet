%Visualize voltage traces - how a spike in one neuron affects neurons it's
%connected to

saveFlag = 0; % 1 to save analysis results

%% ---- import connectivity matrix and voltage traces ----

%Select network to analyze
net_save_path = uigetdir('/Users/hannahgermaine/Documents/PhD/','Select Network Save Folder'); %Have user input where they'd like the output stored
slashes = find(net_save_path == '/');
save_path = net_save_path(1:slashes(end));
scriptFolder = '/voltageVis'; % sub-folder so save analysis results to

if saveFlag & ~isfolder([net_save_path, scriptFolder])
    mkdir([net_save_path, scriptFolder]);
end

%Load data
load(strcat(net_save_path,'/V_m_var.mat')) %want network_var(i).V_m for initialization i
load(strcat(net_save_path,'/network.mat')) %want network(1).conns
load(strcat(save_path,'/parameters.mat')) %want relevant simulation parameters

%Save relevant variables and clear remainder
n = parameters.n; %number of neurons in network
inits = parameters.test_val_max; %number of initializations
V_th = parameters.V_th; %firing threshold
dt = parameters.dt; %timestep (s)
clear parameters
V_m = struct;
for i = 1:inits
    V_m(i).V_m = V_m_var(i).V_m;
end
clear i network_var
conns = network.conns;
clear network

%% ---- select random subset of neurons to look at ----
n_p = 5; %number of primary neurons to look at
n_m = 5; %number of connected neurons to look at

neur_ind = randi(n,[1,n_p]); %which indices of neurons to look at
conn_ind = zeros(n_p,n_m); %which connected neuron indices to look at
for i = 1:n_p
    c_i = find(conns(i,:)); %find indices of connected neurons
    c_i_perm = c_i(randperm(length(c_i))); %randomly permute indices of connected neurons
    if length(c_i) > n_m
        conn_ind(i,:) = c_i_perm(1:n_m); %store only first n_m connected neurons
    else %when there are less connected neurons than n_m
        conn_ind(i,1:length(c_i)) = c_i_perm(1:length(c_i));
    end
end
clear i c_i c_i_perm

%% ---- plot membrane potential deflections ----
%bin = 5*10^(-3); %how long of a deflection to look at (s)
bin_t = 10; %ceil(bin/dt); %how many timesteps of deflection to look at
for i = 1:n_p
    %search for an initialization where neuron i spikes
    for j = 1:inits
        [~,spike_times] = find(V_m(j).V_m >= V_th);
        if ~isempty(spike_times)
            break
        end
    end %j = initialization, spike_times = time indices of spikes
    
    %store the voltage traces of the connected neuron subset
    n_spikes = length(spike_times);
    sub_V_m = zeros(n_spikes,n_m,bin_t+1);
    for k = 1:n_spikes
        spike_t = spike_times(k);
        for l = 1:n_m
            c_i = conn_ind(i,l); %connected neuron index
            if c_i ~= 0 %to account for times of less connected neurons
                sub_V_m(k,l,:) = V_m(j).V_m(c_i,spike_t:spike_t + bin_t);
            end
        end    
    end
    
    %visualize the connected neuron traces as a figure
    f = figure;
    subplot_n = sqrt(n_m);
    if subplot_n == round(subplot_n) %square
        subplot_x = subplot_n;
        subplot_y = subplot_n;
    else
        subplot_x = floor(subplot_n);
        subplot_y = ceil(subplot_n);
    end
    axes = []; %sub_V_m = zeros(n_spikes,n_m,bin_t);
    for l = 1:n_m
        c_i = conn_ind(i,l);
        ax = subplot(subplot_y, subplot_x, l);
        axes = [axes, ax];
        avg_deflection = mean(squeeze(sub_V_m(:,l,:)))*10^3;
        std_deflection = std(squeeze(sub_V_m(:,l,:)))*10^3;
        std1 = avg_deflection + std_deflection;
        std2 = avg_deflection - std_deflection;
        inBetween = [std1, fliplr(std2)];
%         fill([0:bin_t, fliplr(0:bin_t)], inBetween, 'b','FaceAlpha',0.5);
%         hold on
        plot(0:bin_t,avg_deflection,'b')
        xlabel('Timesteps Following Spike')
        ylabel('Membrane Potential (mV)')
        title(strcat('Connected Neuron Index ',string(c_i)))
    end
    %linkaxes(axes)
    sgtitle(strcat('Average Deflections From Spiking of Neuron ',string(neur_ind(i))))
    if saveFlag
        savefig(f,strcat(net_save_path, scriptFolder, '/deflections_from_n',string(neur_ind(i)),'.fig'))
        saveas(f,strcat(net_save_path, scriptFolder, '/deflections_from_n',string(neur_ind(i)),'.jpg'))
        close(f)
    end
end
