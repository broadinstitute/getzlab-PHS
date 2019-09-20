function [L P] = process_results(F)

%takes input of results, outputs loci and p-values in concatenated format

demand_field(F, 'file');

L = cell(slength(F), 1);
P = cell(slength(F), 1);
n_extra = NaN(slength(F), 1);

C = load_struct('ref/1536_LuT.txt');
C = makeapn(C);
C = sparse(C.c512, C.nbidx, C.ch1536);

parfor i = 1:slength(F),
  X = load(F.file{i});

  if ~isfield(X.P, 'pp_recurrence_floor'),
    X.P.pp_recurrence_floor = 0;
  end

  %save mutation struct
  L{i} = reorder_struct(rmfields(X.Mu, {'ref_idx' 'gene_idx'}), ~isnan(X.Mu.p(:, 1)));

  %annotate each locus with Bayes factor
  L{i}.bf = (X.logmarg_lik - X.logmarg_lik_unif)*ones(slength(L{i}), 1);

  %add ch96 annotation
  L{i}.ch96 = F.ch96(i)*ones(slength(L{i}), 1);

  %add ch1536 annotation
  nb = mod(ceil(F.ch96(i)/16) - 1, 3) + 1;
  L{i}.ch1536 = full(C(L{i}.c512, nb));

  %save all p-values
  P{i} = X.pvalues;

  n_extra(i) = nnz(~idx);
end

L = concat_structs(L);
P = cat(1, P{:});
n_extra = sum(n_extra);

L.p(L.p <= 1e-16) = 1e-16;

L.q = NaN(size(L.p));
for i = 1:3,
  L.q(:, i) = fdr_jh(L.p(:, i), n_extra);
end
