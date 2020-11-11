% adopted from:
 %M.J. Beal's implementation of Variational Bayesian Mixture of Factor Analysers:
%http://www.cse.buffalo.edu/faculty/mbeal/software.html
%Ghahramani, Z. and Beal, M.J. (2000)
%Variational Inference for Bayesian Mixtures of Factor Analysers
%In Advances in Neural Information Processing Systems 12:449-455, eds. S. A. Solla, T.K. Leen, K, MIT Press, 2000.

% Modified/adpted by Jalil Taghia

s = size(Lm,2);

if s>1
      temp_mcl = cat(3,Lm{:});
      temp_mcl = squeeze(temp_mcl(:,1,:));
      temp_Lcov = cat(4,Lcov{:,:,:});
      
      mean_mcl = mean(temp_mcl,2);
      nu_mcl = s./( sum(squeeze(temp_Lcov(1,1,:,:)),2) ...
            + sum(temp_mcl.^2,2) ...
            - 2*mean_mcl.*sum(temp_mcl,2) ...
            + s*mean_mcl.^2 );
end
