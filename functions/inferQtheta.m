% adopted from:
 %M.J. Beal's implementation of Variational Bayesian Mixture of Factor Analysers:
%http://www.cse.buffalo.edu/faculty/mbeal/software.html
%Ghahramani, Z. and Beal, M.J. (2000)
%Variational Inference for Bayesian Mixtures of Factor Analysers
%In Advances in Neural Information Processing Systems 12:449-455, eds. S. A. Solla, T.K. Leen, K, MIT Press, 2000.

% Modified/adpted by Jalil Taghia

nSubjs = length(Ycell);
wa = zeros(nStates,nStates); wpi = zeros(1,nStates);
for ns = 1:nSubjs
      phthtpgV1T1 = Qnss{ns}; phtgV1T1 = QnsCell{ns};
      wa = wa + sum(phthtpgV1T1,3);
      wpi = wpi + phtgV1T1(1,:);
end
Wa = wa + repmat(ua,[nStates 1]);
Wpi = wpi + upi;

stran = exp(  psi(Wa) - repmat( psi(sum(Wa,2)) ,[1 nStates])  );
sprior = exp(  psi(Wpi) - psi(sum(Wpi,2))  );