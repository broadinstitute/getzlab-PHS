function [eff_exp_dir eff_obs_dir] = run_effect_counts_dir_mult_MCMC(O, terr, P)

if ~exist('P', 'var'), P = []; end
P = impose_default_value(P, 'a0', 10);
P = impose_default_value(P, 'niter', 5e4);
P = impose_default_value(P, 'collapse_512_to_32', true);

demand_fields(O, {'nb' 'any'});

%collapse to 32 contexts
if P.collapse_512_to_32,
  O_orig = O;
  terr_orig = terr;

  O = [];
  O.nb = NaN(32, 3);
  O.ef = NaN(32, 3);
  O.any = false(32, 1);

  terr = NaN(32, 3, 3);

  for x = [1:16:512; 16:16:512; 1:32],
    i = x(1); j = x(2); k = x(3);

    O.nb(k, :) = nansum(O_orig.nb(i:j, :));
    O.ef(k, :) = nansum(O_orig.ef(i:j, :));
    O.any(k) = any(O_orig.any(i:j));

    terr(k, :, :) = sum(terr_orig(i:j, :, :));
  end
end

F = bsxfun(@rdivide, terr, sum(terr, 2));
F = F(O.any, :, :);
Nc = O.nb(O.any, :);

ne = cell(size(F, 1), 1);
parfor i = 1:size(F, 1),
  ne{i} = effect_counts_dir_mult_MCMC(squeeze(F(i, :, :)), Nc(i, :), P);
end
ne = cat(3, ne{:});
ne_tot = sum(ne, 3);

%distribution on expected effect
eff_exp_dir = NaN(size(ne_tot, 2) - 500, 3);
for i = 1:size(eff_exp_dir, 1),
  eff_exp_dir(i, :) = dirrnd(ne_tot(:, i + 500));
end

%distribution on observed effect
ef_obs = sum(O.ef(O.any, :));
eff_obs_dir = dirrnd(ef_obs, size(eff_exp_dir, 1));
