% adopted from:
 %M.J. Beal's implementation of Variational Bayesian Mixture of Factor Analysers:
%http://www.cse.buffalo.edu/faculty/mbeal/software.html
%Ghahramani, Z. and Beal, M.J. (2000)
%Variational Inference for Bayesian Mixtures of Factor Analysers
%In Advances in Neural Information Processing Systems 12:449-455, eds. S. A. Solla, T.K. Leen, K, MIT Press, 2000.

% Modified/adpted by Jalil Taghia

function [res] = kldirichlet(vecP,vecQ)

alphaP = sum(vecP,2);
alphaQ = sum(vecQ,2);

res = gammaln(alphaP)-gammaln(alphaQ) ...
    - sum(gammaln(vecP)-gammaln(vecQ),2) ...
    + sum( (vecP-vecQ).*(digamma(vecP)-digamma(alphaP)) ,2);