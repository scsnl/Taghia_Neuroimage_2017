function [data,hiddenSources,actual_states, dataCov] = simdata_hmm2(U, transProb, pio, no_subjs, mus, Covs, len,  noise)
if nargin<8
      noise =.01;
end
%Simulate HMM data
for subj = 1:no_subjs
    S = mc_sample(pio, transProb, len);
    y = zeros(size(U,1),len);
    hx = zeros(size(U,2),len);
    states = unique(S);
    for s = 1:length(states)
        ix = find(S == states(s));
        tempx = mvnrnd(mus(states(s),:),Covs(:,:,states(s)),length(ix))';
        temp2 = U(:,:,s) * tempx +repmat(noise.*randn(size(U,1),1),1,size(tempx,2));
        y(:,ix) = temp2;
        dataCov(:,:,s) = cov(temp2');
        hx(:,ix) = tempx;
    end
    data{subj} = y;
    hiddenSources{subj} = hx;
    actual_states{subj} = S;
end