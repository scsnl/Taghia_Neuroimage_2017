% test script

% add the folder ot your path
addpath(genpath('BSFA'));

%%%generate data%%%
nSubjs = 4;
nStates = 4;
nSamps = 200; 
dim = 10;
ldim = dim/2;
[data, hiddenSources, trueStates]  = generate_random_data(nSamps, nSubjs, dim, nStates, ldim);

%%% preprocessing of data to zero mean and unit variance%%%
for subj =1:nSubjs
      data{subj} = zscore(data{subj}')';
end

%%%leanring a group level model%%%
max_nstates = 10; % usually larger than what you expect
max_ldim = dim-1; % in a fully noninformative set-up
% var_percentile = 0.9; 
% max_ldim = bound_on_model_complexity(data, var_percentile);  for high dimensional data
opt.n_iter = 10; % usually less than 200 is enough 
opt.tol = 1e-3; 
opt.noise = 1;
group_model = BayesianSwitchingFactorAnalysis(data, max_nstates, max_ldim, opt);

%%%learning subject-level state-specific stats from trained group model%%%
% if you need subject-level stats use the below function from group model
% [est_covariance, est_mean] = computeSubjectLevelStatsFromGroupLevelModel(data, group_model);

