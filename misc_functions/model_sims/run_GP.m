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

  if max(X.M) == 1, continue; end

  C = [X.C(:, 1:4) full(sparse(1:length(X.C), X.C(:, 5), 1))];
  r = nbreg(C, X.M);

  full_samps{j} = nbinrnd(1/r.alpha, repmat(1./(1 + r.alpha*exp(C*r.b)), 1, 100));
end
pp.delete

%concatenate and save
full_samps = sparse(cat(1, full_samps{:}));

save(sprintf('figures/model_sims/samps_GP_%d.mat', s), 'full_samps', '-v7.3')
