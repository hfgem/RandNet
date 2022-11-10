function plot_randnet_results(parameters, network, V_m, G_in, network_spike_sequences, ithTest, net_save_path)
% 
% Plots basic simulation results from randnet.m
% 

% Formatting options for plotSpikeRaster: uncomment one line below
% MarkerFormat.MarkerSize = 6.7; %MarkerFormat.Marker = '.';
% MarkerFormat.LineWidth = 1.5; MarkerFormat.Marker = '|';
MarkerFormat = struct; % use this line to choose plotSpikeRaster defaults

% Find spikes
spikes_V_m = V_m >= parameters.V_th;
[spikes_x,spikes_t] = find(spikes_V_m);
max_time = max(spikes_t);
spiking_neurons = unique(spikes_x, 'stable');

%% If there are any detected events, plot sorted rasters
if isfield(network_spike_sequences(ithTest), 'events')
    
    events = network_spike_sequences(ithTest).events;

    %Visualize re-ordered spike sequences
    if ~isempty(events)
        
        f = figure;
        num_events = size(events, 1);
        for e_i = 1:num_events

            spike_ranks = network_spike_sequences(ithTest).ranks_vec(:,e_i);
            if numel(spike_ranks)== parameters.n_E % ranks for only E cells
                [~, Ie] = sort(spike_ranks);
                eventSpikes = [spikes_V_m(network.E_indices(Ie),events(e_i,1):events(e_i,2)); ...
                    spikes_V_m(network.I_indices,events(e_i,1):events(e_i,2))];
            elseif numel(spike_ranks)==parameters.n              % ranks for all cells
                [~, Ie] = sort(spike_ranks(network.E_indices));
                [~, Ii] = sort(spike_ranks(network.I_indices));
                eventSpikes = [spikes_V_m(network.E_indices(Ie),events(e_i,1):events(e_i,2)); ...
                    spikes_V_m(network.I_indices(Ii),events(e_i,1):events(e_i,2))];
            else
                error('Unkonwn number of cells in spike_ranks')
            end
            subplot(1,num_events,e_i)
            plotSpikeRaster( eventSpikes, 'TimePerBin', parameters.dt, 'PlotType', 'scatter', 'MarkerFormat', MarkerFormat);
            xlabel('Time (s)','FontSize',16)
            ylabel('Reordered Neuron Number','FontSize',16)
            title(strcat('Event #',string(e_i)))

        end
        sgtitle('Spiking Behavior','FontSize',16)

        %linkaxes(axes)
        if parameters.saveFlag
            savefig(f,strcat(net_save_path,'/','_',string(ithTest),'firing_sequence.fig'))
            saveas(f,strcat(net_save_path,'/', '_',string(ithTest),'firing_sequence.jpg'))
            close(f)
        end
        clear e_i spike_order reordered_spikes event_length s_i ax axes xt xtlbl
    end      
end


%% Plot example Gin and Vm
subSetToPlot = 1:2;
t = [0:parameters.dt:parameters.t_max];
figure; plot(t, V_m(subSetToPlot,:)); ylabel('Vm (V)'); xlabel('Time (s)'); 
figure; plot(t, G_in(subSetToPlot,:)); ylabel('G in (S)'); xlabel('Time (s)'); 


%% Plot population raster
figure; hold on
if exist('events')==1 
    for i = 1:size(events, 1)
        fill([t(events(i,:)), fliplr(t(events(i,:)))], [0, 0, size(spikes_V_m, 1), size(spikes_V_m, 1)], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
    end
end
plotSpikeRaster( [spikes_V_m(network.E_indices,:);  spikes_V_m(network.I_indices,:)], ...
    'TimePerBin', parameters.dt, 'PlotType', 'scatter', 'MarkerFormat', MarkerFormat); % 
ylabel('Cell'); xlabel('Time (s)'); 


%% Plot population firing rate
%{
movmeanWindow = (1/parameters.dt) * 0.05;
meanPopRate = movmean(mean(spikes_V_m, 1)/parameters.dt, movmeanWindow);
figure; hold on
if exist('events')==1 
    for i = 1:size(events, 1)
        fill([t(events(i,:)), fliplr(t(events(i,:)))], [0, 0, ceil(max(meanPopRate)), ceil(max(meanPopRate))], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
    end
end
plot(t, meanPopRate)
ylabel('Population mean rate (Hz)'); xlabel('Time (s)'); 
yline(mean(meanPopRate), 'g')
yline(mean(meanPopRate)+ std(meanPopRate))
yline(mean(meanPopRate)+ 2*std(meanPopRate))
%}

%% Plot E and I cell population firing rate

meanEPopRate = smoothdata(mean(spikes_V_m(network.E_indices,:), 1)/parameters.dt, 'gaussian', parameters.PBE_window);

PBEthresh = max(mean(meanEPopRate)+(parameters.PBE_zscore*std(meanEPopRate)), parameters.PBE_min_Hz); % threshold rate for PBE detection


figure; hold on
if exist('events')==1 
    for i = 1:size(events, 1)
        fill([t(events(i,:)), fliplr(t(events(i,:)))], [0, 0, ceil(max(meanEPopRate)), ceil(max(meanEPopRate))], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
    end
end

% meanIPopRate = movmean(mean(spikes_V_m(network.I_indices,:), 1)/parameters.dt, movmeanWindow);
meanIPopRate = smoothdata(mean(spikes_V_m(network.I_indices,:), 1)/parameters.dt, 'gaussian', parameters.PBE_window);
yyaxis right; 
plot(t, meanIPopRate, ':')
ylabel('I cell mean rate (Hz)'); 

yyaxis left; 
plot(t, meanEPopRate)
ylabel('E cell mean rate (Hz)'); xlabel('Time (s)'); 
yline(mean(meanEPopRate), 'g')
yline(mean(meanEPopRate)+ std(meanEPopRate))
yline(mean(meanEPopRate)+ 2*std(meanEPopRate))
yline(PBEthresh, 'r');



%% Main figure

sortEcells = 2;

if exist('myPlotSettings'); myPlotSettings(8.5, 5.5); end

% 
% myPlotSettings(8.5, 4.5, 2, 14, [], [], 2) % ppt format
% MarkerFormat.MarkerSize = 6;
% myPlotSettings(10, 6, 3, 24, [], [], 3) % SfN-poster format
% MarkerFormat.MarkerSize = 6; % SfN-poster format


figure; 
if ~isempty(events)==1 
    sgtitle({['Frac partic: ', regexprep(num2str( round( mean(~isnan(network_spike_sequences.ranks_vec), 1),  2)),'\s+',',') ], ...
        ['Event dur. (ms): ', regexprep(num2str( round( network_spike_sequences.event_lengths'*1000,  0)),'\s+',',') ] } )  % num2str( mean(~isnan(network_spike_sequences.ranks_vec), 1))
    
    spike_ranks1 = network_spike_sequences.ranks_vec(:,1);
    [~, Ie1] = sort(spike_ranks1);

    
    meanRelRank = nanmean(network_spike_sequences.ranks_vec./parameters.n_E, 2);
    [~, IeRel] = sort(meanRelRank);

else
    sgtitle('No events detected')

    Ie1 = ones(size(network.E_indices));
end

if sortEcells==1 & ~isempty(events)
    reordered_spikes = [spikes_V_m(network.E_indices(Ie1),:); spikes_V_m(network.I_indices,:)];
elseif sortEcells==2 & ~isempty(events) 
    reordered_spikes = [spikes_V_m(network.E_indices(IeRel),:); spikes_V_m(network.I_indices,:)];
else
    reordered_spikes = [spikes_V_m(network.E_indices,:); spikes_V_m(network.I_indices,:)];
end

% Plot Raster
ax1 = subplot(2,1,1); hold on
if exist('events')==1 
    for i = 1:size(events, 1)
        fill([t(events(i,:)), fliplr(t(events(i,:)))], [0, 0, size(spikes_V_m, 1), size(spikes_V_m, 1)], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
    end
end
plotSpikeRaster( reordered_spikes, ...
    'TimePerBin', parameters.dt, 'PlotType', 'scatter', 'MarkerFormat', MarkerFormat); % 
ylabel('Cell');

% Plot population firing rate
ax2 = subplot(2,1,2); hold on
meanEPopRate = smoothdata(mean(spikes_V_m(network.E_indices,:), 1)/parameters.dt, 'gaussian', parameters.PBE_window);
PBEthresh = max(mean(meanEPopRate)+(parameters.PBE_zscore*std(meanEPopRate)), parameters.PBE_min_Hz); % threshold rate for PBE detection
if exist('events')==1 
    for i = 1:size(events, 1)
        fill([t(events(i,:)), fliplr(t(events(i,:)))], [0, 0, ceil(max(meanEPopRate)), ceil(max(meanEPopRate))], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
    end
end
plot(t, meanEPopRate)
ylabel({'Population rate', '(Hz)'}); xlabel('Time (s)'); 

yline(mean(meanEPopRate), 'g')
yline(mean(meanEPopRate)+ std(meanEPopRate))
yline(mean(meanEPopRate)+ 1.5*std(meanEPopRate))
yline(PBEthresh, 'r');

linkaxes([ax1, ax2], 'x')

% For SfN poster
% xlim([3.5, 5.5]); sgtitle '';  ylim([0, 10])
% ylim([0, 375]); set(gca,'XColor','none')


if exist('myPlotSettings'); myPlotSettings; end


end