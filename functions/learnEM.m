inferQtheta;
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