% adopted from:
 %M.J. Beal's implementation of Variational Bayesian Mixture of Factor Analysers:
%http://www.cse.buffalo.edu/faculty/mbeal/software.html
%Ghahramani, Z. and Beal, M.J. (2000)
%Variational Inference for Bayesian Mixtures of Factor Analysers
%In Advances in Neural Information Processing Systems 12:449-455, eds. S. A. Solla, T.K. Leen, K, MIT Press, 2000.

% Modified/adpted by Jalil Taghia

function [kl] = klgamma(pa,pb,qa,qb)

n = max([size(pb,2) size(pa,2)]);

if size(pa,2) == 1, pa = pa*ones(1,n); end
if size(pb,2) == 1, pb = pb*ones(1,n); end
qa = qa*ones(1,n); qb = qb*ones(1,n);

kl = sum( pa.*log(pb)-gammaln(pa) ...
      -qa.*log(qb)+gammaln(qa) ...
      +(pa-qa).*(digamma(pa)-log(pb)) ...
      -(pb-qb).*pa./pb ,2);