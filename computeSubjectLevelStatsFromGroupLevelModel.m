function [est_covariance, est_mean] = computeSubjectLevelStatsFromGroupLevelModel(data, group_model, opt)
%Inputs: 
%	data 				       : {data_subj1, data_subj2,...,data_subjS}, data_subj1 is a D-by-N matrix where D is dimension 
%     a trained BSFA model
%	opt.n_iter 			    : number of iteration of the variational inference
%	opt.tol 			     
%Outputs:	
%	estimated_mean 		  : estimated mean of each state  for each subject
%	estimated_covariance  : estimated covariance of each state for each subject

if nargin<3
	opt.n_iter = 50;
      opt.tol = 1e-3;
end

nSubjs = length(data);
maxLocalDim = size(group_model.net.params.Lm{1}, 2)-1;
for subj = 1:nSubjs
      display(['subject :', num2str(subj)]);
      net_subj = bsfa_z(data(subj), group_model.net.hidden.QnsCell{subj}, maxLocalDim, opt.n_iter,opt.tol , group_model.net.noise);
      est_covariance{subj} = getCovariance(net_subj);
      est_mean{subj} = getMean(net_subj);
end
