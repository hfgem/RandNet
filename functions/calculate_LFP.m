function LFP = calculate_LFP(I_syn, network)
%ABOUT: This function takes in the synaptic current information along with
%which synapses are excitatory and inhibitory and calculates a proxy LFP.
%
%INPUTS:
%   I_syn = matrix of [n x t_steps] with the synaptic current for each
%       neuron (n) across time (t_steps).
%   network = structure file which contains information on which neurons
%   are excitatory and inhibitory.
%
%OUTPUTS:
%   LFP = vector of the estimated LFP

%First find which neurons are post-synaptically pyramidal

pyr_neur = network.E_indices;


%Then sum all incoming excitatory and inhibitory currents ('post-synaptic
%currents') with the following equation:
%   LFP(t) = Norm[sum_pyr(AMPA(t - 6 ms)) - 1.65*sum_pyr(GABA(t))]



%

end