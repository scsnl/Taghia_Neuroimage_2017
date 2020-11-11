function [fractional_occupancy, mean_life]  = compute_occupancy_and_mean_life(temporal_evolution_of_states, max_nstates)

n_subjs = length(temporal_evolution_of_states);
K = max_nstates;
for subj=1:n_subjs
      sc =zeros(1,K);
      for j=1:K
            sc(j) =  sum(temporal_evolution_of_states{subj} ==j);
      end
      counts_post = sc/sum(sc);
      [percent_dominant(subj,:), dominant_states(subj,:)] = sort(counts_post,'descend');
      [fractional_occupancy(subj,:), mean_life(subj,:),~] = summary_stats_fast(temporal_evolution_of_states{subj}, dominant_states(subj,:));
end