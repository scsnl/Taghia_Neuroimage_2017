% adopted from:
 %M.J. Beal's implementation of Variational Bayesian Mixture of Factor Analysers:
%http://www.cse.buffalo.edu/faculty/mbeal/software.html
%Ghahramani, Z. and Beal, M.J. (2000)
%Variational Inference for Bayesian Mixtures of Factor Analysers
%In Advances in Neural Information Processing Systems 12:449-455, eds. S. A. Solla, T.K. Leen, K, MIT Press, 2000.

% Modified/adpted by Jalil Taghia

function [marginalLikelihood, logLikelihood ] = marginalLikelihoodFromGivenNetForNewData(Y, net)


if isempty(net)
error('net obejct is invalid')
end



Ycell = Y;
Y = cell2mat(Y);
psii = net.hparams.psii;
Xm = net.hidden.Xm;
Lm = net.params.Lm;
Xcov = net.hidden.Xcov;
Lcov = net.params.Lcov;
n = size(Y,2);
p = size(Y,1);
s = size(Lm,2);
stran = net.params.stran;
sprior = net.params.sprior;
Wa = net.params.Wa;
Wpi = net.params.Wpi;
Qns = net.hidden.Qns;
a = net.hidden.a;
b = net.hidden.b;
pa = net.hparams.pa;
pb = net.hparams.pb;
mean_mcl = net.hparams.mcl(:,1);
nu_mcl = net.hparams.mcl(:,2);
Fhist = [];
iter = 1;

inferQX;
computeLogOutProbs;

[loglik,QnsCell,Qnss] = vbhmmEstep(Ycell,stran', sprior,logOutProbs);
Qns = [];
for ns=1:nSubjs
    Qns = [Qns; QnsCell{ns}];
end
computeLowerBoundOnPredictiveDensity; 
marginalLikelihood = F;
% logLikelihood = loglik;
logLikelihood = marginalLikelihood;


