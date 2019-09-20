function run_LNP(s)

if isdeployed,
  s = str2double(s);
end

F = [];
F.file = direc('pois_reg/output/*.mat');
F = parsein(F, 'file', '.*output/(\d+)-(.*).mat', {'ch96' 'context'});
F.ch96 = str2double(F.ch96);
F = reorder_struct(F, ~isnan(F.ch96));
F = sort_struct(F, 'ch96');

full_samps = cell(slength(F), 1);

pp = parpool(16);
parfor j = 1:slength(F),
  X = load(F.file{j}); X = X.X;

  %get regression coefficients
  C = [X.C(:, 1:4) full(sparse(1:length(X.C), X.C(:, 5), 1))];
  [BetaP, dBeta, stats] = glmfit(C, X.M, 'poisson', 'constant', 'off');

  full_samps{j} = poissrnd(repmat(exp(C*BetaP), 1, 100));
end
pp.delete

%concatenate and save
full_samps = sparse(cat(1, full_samps{:}));

save(sprintf('figures/model_sims/samps_UP_%d.mat', s), 'full_samps', '-v7.3')
