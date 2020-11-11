function S = mc_sample(prior, trans, len, numex)
if nargin==3
      numex = 1;
end
S = zeros(numex,len);
for i=1:numex
      S(i, 1) = sampleDiscrete(prior);
      for t=2:len
            S(i, t) = sampleDiscrete(trans(S(i,t-1),:));
      end
end

end
function M = sampleDiscrete(prob, r, c)
n = length(prob);
if nargin == 1
      r = 1; c = 1;
elseif nargin == 2
      c = r;
end
R = rand(r, c);
M = ones(r, c);
cumprob = cumsum(prob(:));
if n < r*c
      for i = 1:n-1
            M = M + (R > cumprob(i));
      end
else
      cumprob2 = cumprob(1:end-1);
      for i=1:r
            for j=1:c
                  M(i,j) = sum(R(i,j) > cumprob2)+1;
            end
      end
end

end
