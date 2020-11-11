% adopted from:
 %M.J. Beal's implementation of Variational Bayesian Mixture of Factor Analysers:
%http://www.cse.buffalo.edu/faculty/mbeal/software.html
%Ghahramani, Z. and Beal, M.J. (2000)
%Variational Inference for Bayesian Mixtures of Factor Analysers
%In Advances in Neural Information Processing Systems 12:449-455, eds. S. A. Solla, T.K. Leen, K, MIT Press, 2000.

% Modified/adpted by Jalil Taghia
s = size(Lm,2);
a = pa + .5*size(Lm{1},1);
for tt = 1:s
  kt = size(Lm{tt},2);
  b{tt} = pb*ones(1,kt-1) + .5*( diag(sum(Lcov{tt}(2:end,2:end,:),3))' + sum(Lm{tt}(:,2:end).^2,1) );
end


