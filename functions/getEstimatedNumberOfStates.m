function estimated_number_of_states = getEstimatedNumberOfStates(temporal_evolution_of_states)

n_subjs = length(temporal_evolution_of_states);
for subj = 1:n_subjs
      estimated_number_of_states{subj} = unique(temporal_evolution_of_states{subj});
end