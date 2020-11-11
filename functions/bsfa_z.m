% adopted from:
 %M.J. Beal's implementation of Variational Bayesian Mixture of Factor Analysers:
%http://www.cse.buffalo.edu/faculty/mbeal/software.html
%Ghahramani, Z. and Beal, M.J. (2000)
%Variational Inference for Bayesian Mixtures of Factor Analysers
%In Advances in Neural Information Processing Systems 12:449-455, eds. S. A. Solla, T.K. Leen, K, MIT Press, 2000.

% Modified/adpted by Jalil Taghia

function net = bsfa_z(Y, Qns, maxdim, nIter,tol, pcaflag)
Ycell = Y;
Y = cell2mat(Y);

nStates = size(Qns, 2);
[p, n] = size(Y);
k = maxdim + 1;
net = struct( ...
      'model','variational SFA_GivenLatent_variables', ...
      'hparams',struct('mcl',[],'psii',[],'pa',[],'pb',[],...
      'alpha',[],'alphaa',[],'alphapi',[]), ...
      'params',struct('Lm',[],'Lcov',[],'u',[],'stran',[],'sprior',[]), ...
      'hidden',struct('Xm',[],'Xcov',[],'Qns',[],'Qnss',[],'a',[],'b',[]), ...
      'Fhist',[]);
pa = 1; pb = 1; alpha=1; alphaa = 1; alphapi = 1; psii = ones(p,1); u = 1;
for ss=1:nStates
      Lm{ss} = randn(p,k); Lcov{ss} = repmat(eye(k),[1 1 p]);
      mean_mcl = mean(Y,2); nu_mcl = (1./std(Y,0,2)).^2;
      Xcov{ss} = zeros(k,k,n);
end
psimin = 1e-5;

s = nStates;
F = -inf;
Fhist = [];
nSubjs = length(Ycell);
% learning for the first iteration
inferQnu;
inferQX;
inferQL;
inferpsii2;
infermcl;
computeLogOutProbs;
iter = 1;
ll = [];
F = -inf;
dF = inf;
ddF = inf;
criteria = 1;
reverseStr = '';
while  criteria
      learnEM_z;
      computeLowerBound_z;
      if s~=1
            ll = [ll F];
      end
      if iter==1
            dF = Inf;
      else
            temp = abs(abs(ll(end))-abs(ll(end-1)));
            dF = [dF temp(end)];
            ddF = abs(abs(dF(end))-abs(dF(end-1)));
      end
      iter = iter+1;
      if isinf(abs(dF(end))) || isnan(abs(dF(end))) || isinf(abs(ddF)) || isnan(abs(ddF))
            criteria = iter<nIter ||  (abs(dF(end))>tol  &&  abs(ddF)>tol);
      else
            criteria = iter<nIter &&  (abs(dF(end))>tol  &&  abs(ddF)>tol);
      end
      msg = sprintf('Processed %d/%d', iter, nIter);
      fprintf([reverseStr, msg]);
      reverseStr = repmat(sprintf('\b'), 1, length(msg));
end
display('  ')
net.hparams.mcl = [mean_mcl nu_mcl];
net.hparams.psii = psii;
net.hparams.alpha = alpha;
net.hparams.alpha = alpha;
net.params.Lm = Lm;
net.params.Lcov = Lcov;
net.params.u = u;
net.hidden.Xm = Xm;
net.hidden.Xcov = Xcov;
net.hidden.a = a;
net.hidden.b = b;
net.Fhist =  Fhist;
net.logOutProbs = logOutProbs;



