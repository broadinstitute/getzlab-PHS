function run_LNP(s)

if isdeployed,
  s = str2double(s);
end

F = [];
F.file = direc('models_v3/20171125_a212692/96_regr_hier/PANCAN/output/*.mat');
F = parsein(F, 'file', '.*output/(\d+)-(.*)_3.mat', {'ch96' 'context'});
F.ch96 = str2double(F.ch96);
F = sort_struct(F, 'ch96');

full_samps = cell(slength(F), 1);

pp = parpool(16);
parfor j = 1:slength(F),
  X = load(F.file{j});

  %randomly sample 20 draws from the posterior to average over
  rs = 500 + randsample(2500, 20);

  %sample from posterior
  sig_tots = sqrt(bsxfun(@plus, 1./X.tau(1:16, :), 1./X.tau(end, :)));

  catM = full(sparse(1:size(X.C, 1), X.C(:, 5), 1));

  betaC = X.C(:, 1:4)*X.Beta(:, rs);
  mu = catM*X.mu(:, rs);
  sig = catM*sig_tots(:, rs);

  %draws from posterior predictive (100 total samples per locus)
  samps = NaN([size(betaC) 5]);
  for i = 1:5,
    samps(:, :, i) = poissrnd(exp(normrnd(betaC + mu + X.P.log_exposure, sig)));
  end
  full_samps{j} = reshape(samps, size(betaC, 1), [], 1);
end
pp.delete

%concatenate and save
full_samps = sparse(cat(1, full_samps{:}));

save(sprintf('figures/model_sims/samps_LNP_%d.mat', s), 'full_samps', '-v7.3')
