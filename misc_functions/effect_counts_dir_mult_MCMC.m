function [ne b] = effect_counts_dir_mult_MCMC(F, Nc, P)

if ~exist('P', 'var'), P = []; end
P = impose_default_value(P, 'niter', 1e4);
P = impose_default_value(P, 'a0', 120);

Nc = as_column(Nc); 
N = sum(Nc);

ne = NaN(3, P.niter);
b = NaN(3, P.niter);

%initialize parameters
b(:, 1) = [1 1 1]'/3;
ne(:, 1) = round(F*Nc);

acc_b = NaN(P.niter, 1);

for i = 2:P.niter,
  %M-H sample p(b|-)
  [b(:, i) acc_b(i)] = bsamp(b(:, i - 1), ne(:, i - 1), F, Nc, P.a0);

  %gibbs sample p(ne|-)
  ne(:, i) = mnrnd(N, F*b(:, i));
end

end

function [b acc_b] = bsamp(b, ne, F, Nc, a0)
  %propose from Dirichlet distribution with mode at current value of b
  bP = dirrnd(b*a0 - 3*b + 1)';

  %log posterior ratio
  prat = nansum(ne.*log(F*bP) + Nc.*log(bP)) - nansum(ne.*log(F*b) + Nc.*log(b));

  %log proposal ratio
  aP = (bP*a0 - 3*bP); a = (b*a0 - 3*b);
  qrat = aP'*log(b) + gammaln(sum(aP + 1)) - sum(gammaln(aP + 1)) - ...
         (a'*log(bP) + gammaln(sum(a + 1)) - sum(gammaln(a + 1)));

  %M-H acceptance test
  acc_b = 0;
  if log(rand) < min(0, prat + qrat), b = bP; acc_b = 1; end
end
