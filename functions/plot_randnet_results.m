function plot_randnet_results(parameters, network, V_m, G_in, network_spike_sequences, ithTest, net_save_path)

%Find spike profile
spikes_V_m = V_m >= parameters.V_th;
[spikes_x,spikes_t] = find(spikes_V_m);
max_time = max(spikes_t);
spiking_neurons = unique(spikes_x, 'stable');

if isfield(network_spike_sequences(ithTest), 'events')
    events = network_spike_sequences(ithTest).events;

    %Visualize re-ordered spike sequences
    if any(strcmp('spike_order',fieldnames(network_spike_sequences)))
        f = figure;
        axes = [];
        num_events = size(events, 1)
        for e_i = 1:num_events

            spike_ranks = network_spike_sequences(ithTest).spike_ranks.(strcat('sequence_',string(e_i)));
            if parameters.E_events_only
                [~, Ie] = sort(spike_ranks);
                eventSpikes = [spikes_V_m(network.E_indices(Ie),events(e_i,1):events(e_i,2)); ...
                    spikes_V_m(network.I_indices,events(e_i,1):events(e_i,2))];
            else
                [~, Ie] = sort(spike_ranks(network.E_indices));
                [~, Ii] = sort(spike_ranks(network.I_indices));
                eventSpikes = [spikes_V_m(network.E_indices(Ie),events(e_i,1):events(e_i,2)); ...
                    spikes_V_m(network.I_indices(Ii),events(e_i,1):events(e_i,2))];
            end
            subplot(1,num_events,e_i)
            plotSpikeRaster( eventSpikes, 'TimePerBin', parameters.dt, 'PlotType', 'scatter');
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

t = [0:parameters.dt:parameters.t_max];
figure; plot(t, V_m(1:2,:)); ylabel('Vm (V)'); xlabel('Time (s)'); 
% figure; plot(t, V_m); ylabel('Vm (V)'); xlabel('Time (s)'); 
figure; plot(t, G_in(1:2,:)); ylabel('G in (S)'); xlabel('Time (s)'); 
% figure; plot(t, G_in(1:2,:)); ylabel('G in (S)'); xlabel('Time (s)'); 

figure; hold on
if exist('events')==1 
    for i = 1:size(events, 1)
        fill([t(events(i,:)), fliplr(t(events(i,:)))], [0, 0, size(spikes_V_m, 1), size(spikes_V_m, 1)], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
    end
end
plotSpikeRaster( [spikes_V_m(network.E_indices,:);  spikes_V_m(network.I_indices,:)], 'TimePerBin', parameters.dt, 'PlotType', 'scatter'); 
ylabel('Cell'); xlabel('Time (s)'); 


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

end