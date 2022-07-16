function [event_lengths, event_counts] = pull_avalanches(used_spikes_mat)
    %This function is used to pull out segments of firing from a network to
    %then test for criticality.

    event_lengths = [];
    event_counts = [];
    
    [~,spikes_t] = find(used_spikes_mat);

    if ~isempty(spikes_t)
        last_start = 1;
        last_time = 1;
        spike_counts = 1;
        for t_i = 2:length(spikes_t)
            s_i = spikes_t(t_i);
            if s_i - last_time <= 1 %check if the next spike is in the same or next time bin
                last_time = s_i; %store the updated last time of a spike in the event
                spike_counts = spike_counts + 1;
            else %if the spike is not within an event timebin update parameters to store a new event
                e_len = last_time - last_start;
                event_lengths = [event_lengths, e_len]; %#ok<AGROW> %add the last range of spikes to the events vector
                event_counts = [event_counts, spike_counts]; %#ok<AGROW>
                last_start = s_i;
                last_time = s_i;
                spike_counts = 1;
            end
        end
    end

end