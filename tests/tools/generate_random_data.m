function [data, hiddenSources, trueStates]  = generate_random_data(nSamps, nSubjs, dim, nStates, ldim)

pio = rand(nStates,1); pio = pio./sum(pio);
transProb = random('uni',0.5,0.82,nStates,nStates) + (diag(random('uni',1,2,nStates,1)).*eye(nStates));
transProb = (transProb./repmat(sum(transProb ,2),1,nStates))';
mus = randn(nStates,ldim);
Covs = zeros(ldim,ldim,nStates);
for state=1:nStates
      Covs(:,:,state) = diag(rand(1,ldim));
end
U = zeros(dim,ldim,nStates);
for state=1:nStates
      tempU = randn(dim,ldim);
      U(:,:,state) = tempU - repmat(mean(tempU,2), 1,ldim);
end
[data, hiddenSources,trueStates, ~] = simdata_hmm2(U,transProb,pio,nSubjs,mus,Covs,nSamps);
