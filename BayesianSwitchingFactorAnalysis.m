function model = BayesianSwitchingFactorAnalysis(data, max_nstates, max_ldim, opt)
%Inputs: 
%   data_subj1                                     : {data_subj1, data_subj2,...,data_subjS}, data_subj1 is a D-by-N matrix where D is dimension 
%                                                     and N is the number of samples
%   max_nstates                                    : maximum number of states (set to a value greater than the expected number of states)
%   max_ldim                                       : bound on the latent space dimensionality. default: max_ldim = D-1. max_ldim must be smaller 
%                                                      than dimensionality of data
%   opt.n_iter                                     : number of iteration of the variational inference
%   opt.noise                                      : type of measurement noise, can be either 1 or 0
%                                                    if opt.noise=1, noise variance the same across dimensions
%                                                    if opt.noise=0 (default), noise variance is allowed to vary across dimensions

%Outputs:			
%   model.net                                            : posterior parameters over all model parameters
%   model.estimated_mean 		                         : estimated mean of each state (group-level)
%   model.estimated_covariance                           : estimated covariance of each state (group-level)
%   model.posterior_probabilities                        : posterior probabilities of states across time for each subject q(Z)
%   model.joint_posterior_probabilities	                 : joint posterior probabilities of states across time for each subject q(Z_n-1, Z)
%   model.state_transition_probabilities                 : HMM estimated transition probabilities across states 
%   model.temporal_evolution_of_states                   : estimated temporal evolution of states
%   model.fractional_occupancy_group_wise                : occupancy rate of each state computed in group sense
%   model.mean_lifetime_group_wise                       : mean life time of each state computed in group sense
%   model.id_of_dominant_states_group_wise               : id of dominant states group wise
%   model.fractional_occupancy_subject_wise              : occupancy rate of each state computed in subject sense
%   model.mean_lifetime_subject_wise                     : mean life time of each state computed in subject sense
%   model.id_of_dominant_states_subject_wise             : id of dominant states subject wise
%   model.id_of_remaining_states                         : id of remaining states for each subject
 
if nargin<3
	max_ldim = size(data{1}, 1) - 1;
	opt.n_iter = 100;
      opt.tol = 1e-3;
	opt.noise = 0;
end

if nargin<4
	opt.n_iter = 100;
      opt.tol = 1e-3;
	opt.noise = 0;
end

% modeling 
display('learning...')  
net = bsfa(data, max_nstates, max_ldim, opt.n_iter, opt.tol, opt.noise);

% computing temporal evolution of states
display('computing temporal evolution of states...')  
[~, stateCell] =  estimateStatesByVitterbi(data,net.params,net.logOutProbs);
temporal_evolution_of_states = stateCell;
id_of_remaining_states = getRemainingStateIds(temporal_evolution_of_states);
id_of_dominant_states_group = getDominantStateIdsGroup(temporal_evolution_of_states, max_nstates);
id_of_dominant_states_subject = getDominantStateIdsSubject(temporal_evolution_of_states, max_nstates);

% compute mean life and occupancy rate
display('computing mean life and fractional occupancy...')  
[fractional_occupancy_group, mean_life_group]  = compute_occupancy_and_mean_life_group_wise(temporal_evolution_of_states, max_nstates);
[fractional_occupancy_subj, mean_life_subj]  = compute_occupancy_and_mean_life_subject_wise(temporal_evolution_of_states, max_nstates);

% computing mean and covariance from model
display('computing mean and covariance...')  
estimated_covariance = getCovariance(net);
estimated_mean = getMean(net);

model.net = net;
model.estimated_covariance = estimated_covariance;
model.estimated_mean = estimated_mean;
model.temporal_evolution_of_states = temporal_evolution_of_states;
model.posterior_probabilities = model.net.hidden.Qns;
model.joint_posterior_probabilities = model.net.hidden.Qnss;
model.fractional_occupancy_group_wise = fractional_occupancy_group;
model.mean_lifetime_group_wise = mean_life_group;
model.fractional_occupancy_subject_wise = fractional_occupancy_subj;
model.mean_lifetime_subject_wise = mean_life_subj;
model.state_transition_probabilities = model.net.params.stran;
model.id_of_dominant_states_group_wise = id_of_dominant_states_group;
model.id_of_dominant_states_subject_wise = id_of_dominant_states_subject;
model.id_of_remaining_states = id_of_remaining_states;
