% adopted from:
 %M.J. Beal's implementation of Variational Bayesian Mixture of Factor Analysers:
%http://www.cse.buffalo.edu/faculty/mbeal/software.html
%Ghahramani, Z. and Beal, M.J. (2000)
%Variational Inference for Bayesian Mixtures of Factor Analysers
%In Advances in Neural Information Processing Systems 12:449-455, eds. S. A. Solla, T.K. Leen, K, MIT Press, 2000.

% Modified/adpted by Jalil Taghia

function net = bsfa(Y, nStates, maxdim, nIter,tol,pcaflag,net)
Ycell = Y;
Y = cell2mat(Y);
[p n] = size(Y);

if nargin<2
      error('at least two inputes are expected: data, number of States')
end
if nargin<3,
      k = (p-1) + 1; 
else
      k = maxdim + 1;
end
if nargin<4, nIter = 10; end
if nargin<5, tol = 1e-3; end
if nargin<6, pcaflag = 0; end

if nargin==7 
      Fhist = net.Fhist;
      mean_mcl = net.hparams.mcl(:,1);
      nu_mcl = net.hparams.mcl(:,2);
      psii = net.hparams.psii;
      pa = net.hparams.pa;
      pb = net.hparams.pb;
      alpha = net.hparams.alpha;
      if ~isempty(net.params)
            Lm = net.params.Lm;
            Lcov = net.params.Lcov;
            u = net.params.u;
      else 
            k = (p-1) +1;
            Lm{1} = randn(p,k); Lcov{1} = repmat(eye(k),[1 1 p]); u = 1;
      end
      net.hidden = [];
else 
      net = struct( ...
            'model','variational MFA', ...
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
end
psimin = 1e-5;
s = nStates;
pu = alpha*ones(1,s)/s;
ua = ones(1,s)*(alphaa/s);
upi = ones(1,s)*(alphapi/s);
F = -inf;
Fhist = [];
[Wa,Wpi,QnsCell] = initPoteriors(Ycell,nStates);
nSubjs = length(Ycell);
Qns = [];
for ns=1:nSubjs
      Qns = [Qns; QnsCell{ns}];
end
stran = exp(  psi(Wa) - repmat( psi(sum(Wa,2)) ,[1 s])  );
sprior = exp(  psi(Wpi) - psi(sum(Wpi,2))  );
% learning for the first iteration
inferQnu;
inferQX;
inferQL;
inferpsii2;
infermcl;
computeLogOutProbs;
[loglik,QnsCell,Qnss] = vbhmmEstep(Ycell,stran', sprior,logOutProbs);
Qns = [];
for ns=1:nSubjs
      Qns = [Qns; QnsCell{ns}];
end

iter = 1;
ll = [];
F = -inf;
dF = inf;
ddF = inf;
criteria = 1;
reverseStr = '';
current_weights = -inf*sum(Qns) ;
improvement = inf;

while  criteria
      learnEM;
      computeLowerBound;
      if s~=1
            ll = [ll F+loglik];
      else
            ll = [ll F];
      end
      if iter==1
            dF = Inf;
      else
            temp = abs(abs(ll(end))-abs(ll(end-1)));
            dF = [dF temp(end)];
            ddF = abs(abs(dF(end))-abs(dF(end-1)));
      end
      improvement = sum(abs(sum(Qns) - current_weights));
      iter = iter+1;
      current_weights = sum(Qns);
      if isinf(abs(dF(end))) || isnan(abs(dF(end))) || isinf(abs(ddF)) || isnan(abs(ddF))
            criteria =  (improvement>tol && iter<nIter) ||  (abs(dF(end))>tol  &&  abs(ddF)>tol);
      else
            criteria = iter<nIter &&  (abs(dF(end))>tol  &&  abs(ddF)>tol) &&  improvement>tol;
      end
      msg = sprintf('Processed %d/%d', iter, nIter);
      fprintf([reverseStr, msg]);
      reverseStr = repmat(sprintf('\b'), 1, length(msg));
end
display('  ')
display(['needed_iter= ', mat2str(iter), '  .......  dF= ', mat2str(dF(end)), '  .......  ddF= ', mat2str(ddF(end)), '  .......  improvement= ', mat2str(num2str(improvement))])
display(['wieghts= ', mat2str(sum(Qns))]);

% normalize to account for the epsilon in the computation of HMM
sum_stran = sum(stran,2);
for state=1:nStates
      stran(state, :) = stran(state,:)/sum_stran(state);
end

net.hparams.mcl = [mean_mcl nu_mcl];
net.hparams.psii = psii;
net.hparams.pa = pa;
net.hparams.pb = pb;
net.hparams.alpha = alpha;
net.hparams.alpha = alpha;
net.params.stran = stran;
net.params.sprior = sprior;
net.params.Wa = Wa;
net.params.Wpi = Wpi;
net.params.Lm = Lm;
net.params.Lcov = Lcov;
net.params.u = u;
net.hidden.Xm = Xm;
net.hidden.Xcov = Xcov;
net.hidden.Qns = Qns;
net.hidden.QnsCell = QnsCell;
net.hidden.a = a;
net.hidden.b = b;
net.hidden.Qns = Qns;
net.hidden.Qnss = Qnss;
net.Fhist =  Fhist;
net.logOutProbs = logOutProbs;
net.noise = pcaflag;




