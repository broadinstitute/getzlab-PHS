%
% 1. Log-normal Poisson model {{{

%
%compile

% cd LNP
% mcc -m -d mcc \
% `find ref -maxdepth 1 -type f -exec echo -n "-a {} " \;` \
% `find src -type d -exec echo -n "-I {} " \;` \
% -I ../funcs \
% src/pois_LN_reg_wrapper

%
% generate run files {{{

%
% pan-cancer {{{
clear

mafpath = 'mutation_data/MC3.align75.ICE_PoN-uniqued.M';

basedir = 'LNP_posteriors/pancancer';
ede(basedir)

%base params file
load('LNP_posteriors/params.mat', 'P')

P.covar_dir = {'ref/cod_dnasehs_covars.zerofilled.align75_filtered/v3' ...
               'ref/cod_nuc_covars.align75_filtered/Gm12878' ...
               'ref/cod_mutsig_transformed_covars.align75_filtered/exprmax_transformed' ...
               'ref/cod_mutsig_transformed_covars.align75_filtered/rt_transformed' ...
               'ref/cod_c512_covars.align75_filtered/v1'};
P.compute_posterior_predictive = 1;
P.pp_recurrence_floor = 1;
P.compute_marginal_likelihood = 0;
P.init_from_MAP = 0;
P.tailor_hyperparameters = 1;
P.niter = 3e3;
P.random_effect_cat_covar_index = 'last';
P.skip_figures = 1;
P.log_exposure = log(9019); %original number of patients was 9023, but four patients had mutations
                            %solely in regions we filtered out.
P.outdir = [basedir '/output'];
ede(P.outdir)

%RC LuT
c1024 = dec2base(0:1023, 4) - 48;
c1024 = [(1:1024)' sum(bsxfun(@times, 3 - c1024(:, [1 3 2 5 4]), [4^4 4^3 4^2 4 1]), 2) + 1];

%generate ch96 names
bc = {'A->C' 'A->G' 'A->T' 'C->A' 'C->G' 'C->T'};
B = {'A' 'C' 'G' 'T'};

tmp = dec2base((1:96) - 1, 4) - 48;
name96 = strcat(B(tmp(:, end - 1) + 1), '(', bc(sum(bsxfun(@times, tmp(:, 1:end - 2), [4 1]), 2) + 1), ')', B(tmp(:, end) + 1))';

%enumerate pentamers associated with each trimer+basechange
[c512 nb] = meshgrid(1:512, 1:3);
c512_nb = [c512(:) nb(:)];

Y = [];
Y.c512 = c512_nb(:, 1); 
b = dec2base(Y.c512 - 1, 4, 5) - 48;

Y.nb = c512_nb(:, 2);
Y.ch96 = 1 + [b(:, [2 3]) Y.nb - 1 b(:, 1)]*[1 4 16 48]';

Y = sort_struct(Y, {'ch96'});
[~, ui] = unique(Y.ch96);

%loop over ch96's, generate runfile
f = fopen([basedir '/runs.cmd'], 'w');
for x = [ui [ui(2:end) - 1; slength(Y)]]',
  i = x(1); j = x(2);

  P.output_name = sprintf('%d-%s_%d', Y.ch96(i), name96{Y.ch96(i)}, f);
  paramspath = sprintf('%s/params_%d.mat', basedir, Y.ch96(i));
  save(paramspath, 'P') 

  contexts = strjoin(num2cellstr([Y.c512(i:j); c1024(Y.c512(i:j), 2)]), ',');

  fprintf(f, 'mcc/pois_LN_reg_wrapper %s %s c1025 %d %s\n', mafpath, contexts, Y.nb(i), paramspath);
end
fclose(f);

%}}}

%
% signature subcohorts {{{
clear

%
%index channels of interest
X = generate_signature_subcohort_definitions('mutation_data/MC3.align75.ICE_PoN-uniqued.sig_cohorts/*.M');

%RC LuT
c1024 = dec2base(0:1023, 4) - 48;
c1024 = [(1:1024)' sum(bsxfun(@times, 3 - c1024(:, [1 3 2 5 4]), [4^4 4^3 4^2 4 1]), 2) + 1];

%base params file
load('LNP_posteriors/params.mat', 'P')

covar_blocks = {{'ref/cod_dnasehs_covars.zerofilled.align75_filtered/v3' 'ref/cod_nuc_covars.align75_filtered/Gm12878'} ...
                {'ref/cod_mutsig_transformed_covars.align75_filtered/exprmax_transformed' 'ref/cod_mutsig_transformed_covars.align75_filtered/rt_transformed'}};

P = rmfield(P, 'covar_names');

P.compute_posterior_predictive = 0;
P.compute_marginal_likelihood = 0;
P.init_from_MAP = 0;
P.interval_list = 'ref/target_list.align75_filtered_tier4.txt';
P.tailor_hyperparameters = 1;
P.niter = 3e3;
P.skip_figures = 1;

%generate ch96 names
bc = {'A->C' 'A->G' 'A->T' 'C->A' 'C->G' 'C->T'};
B = {'A' 'C' 'G' 'T'};

tmp = dec2base((1:96) - 1, 4) - 48;
name96 = strcat(B(tmp(:, end - 1) + 1), '(', bc(sum(bsxfun(@times, tmp(:, 1:end - 2), [4 1]), 2) + 1), ')', B(tmp(:, end) + 1))';

%
%generate runfiles (trimers) {{{

%
% runfiles for each signature (NMF) {{{

%overwrite previous definition of X with one tailored to NMF output

X = generate_signature_subcohort_definitions('mutation_data/MC3.align75.ICE_PoN-uniqued.sig_cohorts_NMF/*.M');

%loop over signature subcohorts

for i = 1:slength(X),
  mafpath = X.mafpath{i};
  M = loadM(mafpath);

  basedir = ['LNP_posteriors/signature_subcohorts/' X.names{i}];
  ede(basedir) 
  ede([basedir '/params'])

  %P.log_exposure = log(X.npat(i));
  P.outdir = [basedir '/output/'];
  ede(P.outdir)
  
  %8 runfiles: for each combination of covariate blocks/whether we use pentamers
  fs = NaN(8, 1);
  for f = 1:8,
    fs(f) = fopen(sprintf('%s/runs_%d.cmd', basedir, f - 1), 'w');
  end

  for j = 1:length(X.output{i}),
    Y = [];
    Y.c512 = X.output{i}{j}(:, 1); 
    b = dec2base(Y.c512 - 1, 4, 5) - 48;

    Y.nb = X.output{i}{j}(:, 2);
    Y.ch96 = 1 + [b(:, [2 3]) Y.nb - 1 b(:, 1)]*[1 4 16 48]';

    Y = sort_struct(Y, {'ch96'});
    [~, ui] = unique(Y.ch96);

    for x = [ui [ui(2:end) - 1; slength(Y)]]',
      k = x(1); l = x(2);

      %skip if too few mutations for meaningful fit
      %50 is empirically chosen cutoff.
      if nnz(M.mut.count_nb(M.mut.tier == 4 & ismember(M.mut.c512, Y.c512(k:l)), Y.nb(k)) >= 2) < 50, continue; end

      contexts = strjoin(num2cellstr([Y.c512(k:l); c1024(Y.c512(k:l), 2)]), ',');

      %number of patients for this channel
      P.log_exposure = log(full(M.npat(Y.ch96(k))));

      %write output files
      for f = 0:7,
	fidx = logical(dec2base(mod(f, 4), 2, 2) - 48);
	P.covar_dir = cat(2, covar_blocks{fidx});

	%whether we use context as covariates
	if f <= 3,
	  P.covar_dir = cat(2, P.covar_dir, 'ref/cod_c512_covars.align75_filtered/v1');
	  P.random_effect_cat_covar_index = 'last';
	else
	  if isfield(P, 'random_effect_cat_covar_index'),
	    P = rmfield(P, 'random_effect_cat_covar_index');
	  end
	end

	if isempty(P.covar_dir), P.covar_dir = 'none'; end

	paramspath = sprintf('%s/params/params_%d_%d.mat', basedir, Y.ch96(k), f);
	P.output_name = sprintf('%d-%s_%d', Y.ch96(k), name96{Y.ch96(k)}, f);

	save(paramspath, 'P')

	fprintf(fs(f + 1), 'mcc/pois_LN_reg_wrapper %s %s c1025 %d %s\n', mafpath, contexts, Y.nb(k), paramspath);
      end
    end
  end

  for f = 1:8, fclose(fs(f)); end
end

% }}}

%
% runfiles for UV signature subcohort redux {{{

% while we're at it, why not use the NMF-defined UV subcohort just to be consistent
% this also keeps the output path names unique

X = generate_signature_subcohort_definitions('mutation_data/MC3.align75.ICE_PoN-uniqued.sig_cohorts_NMF/*.M');

i = 7;

mafpath = X.mafpath{i};
M = loadM(mafpath);

basedir = ['LNP_posteriors/signature_subcohorts/' X.names{i}];
ede(basedir) 
ede([basedir '/params'])

P.outdir = [basedir '/output/'];
ede(P.outdir)

%get context512 names
LuT = load_struct('ref/context_1025_categs.txt');
LuT = makeapn(LuT);
LuT = reorder_struct_exclude(LuT, 1025);
LuT.name = cellfun(@(x) [x(6:7) x(1) x(9:10)], LuT.name, 'unif', 0);

%only one regex necessary for UV signature
j = 1;

Y = [];
Y.c512 = X.output{i}{j}(:, 1); 
b = dec2base(Y.c512 - 1, 4, 5) - 48;

Y.nb = X.output{i}{j}(:, 2);
Y.ch96 = 1 + [b(:, [2 3]) Y.nb - 1 b(:, 1)]*[1 4 16 48]';

Y = sort_struct(Y, {'ch96'});
[~, ui] = unique(Y.ch96);

%init the two runfiles' handles
fs = NaN(9, 1);
for f = 8:9,
  fs(f) = fopen(sprintf('%s/runs_%d.cmd', basedir, f), 'w');
end

for x = [ui [ui(2:end) - 1; slength(Y)]]',
  k = x(1); l = x(2);

  contexts = [Y.c512(k:l) c1024(Y.c512(k:l), 2)];

  %number of patients for this channel
  P.log_exposure = log(full(M.npat(Y.ch96(k))));

  %write to both output files (XR + intrinsic covars, intrinsic covars only)
  for f = 8:9,
    if f == 8,
      P.covar_dir = cat(2, covar_blocks{:}, direc('ref/cod_XR_covars.align75_filtered/XPC*')');
    else
      P.covar_dir = cat(2, covar_blocks{:});
    end

    %loop over pentamer contexts
    for q = 1:size(contexts, 1),
      c1024s = contexts(q, :);

      %skip if too few mutations for meaningful fit
      %50 is empirically chosen cutoff.
      if nnz(M.mut.count_nb(M.mut.tier == 4 & ismember(M.mut.c512, c1024s(1)), Y.nb(k)) >= 2) < 50, continue; end

      paramspath = sprintf('%s/params/params_%d_%d.mat', basedir, c1024s(1), f);

      %note that this output nomenclature is inconsistent with everything else, since a separate 
      %output is saved for each pentamer.  since there is only one newbase for UV, we will index 
      %by c512 only.
      %P.output_name = sprintf('%d-%s_%d', Y.ch96(k), name96{Y.ch96(k)}, 8);
      P.output_name = sprintf('%d-%s_%d', c1024s(1), LuT.name{c1024s(1)}, f);

      save(paramspath, 'P')

      fprintf(fs(f), 'mcc/pois_LN_reg_wrapper %s %s c1025 %d %s\n', mafpath, strjoin(num2cellstr(c1024s), ','), Y.nb(k), paramspath);
    end
  end
end

fclose(fs(8));
fclose(fs(9));

% }}}

%
% runfiles for APOBEC hairpin {{{

X = generate_signature_subcohort_definitions('mutation_data/MC3.align75.ICE_PoN-uniqued.sig_cohorts_NMF/*.M');

i = 1;

mafpath = X.mafpath{i};
M = loadM(mafpath);

basedir = ['LNP_posteriors/signature_subcohorts/' X.names{i}];
ede(basedir) 
ede([basedir '/params'])

P.outdir = [basedir '/output/'];
ede(P.outdir)

%get context512 names
LuT = load_struct('ref/context_1025_categs.txt');
LuT = makeapn(LuT);
LuT = reorder_struct_exclude(LuT, 1025);
LuT.name = cellfun(@(x) [x(6:7) x(1) x(9:10)], LuT.name, 'unif', 0);

%only one regex necessary for APOBEC signature
j = 1;

Y = [];
Y.c512 = X.output{i}{j}(:, 1); 
b = dec2base(Y.c512 - 1, 4, 5) - 48;

Y.nb = X.output{i}{j}(:, 2);
Y.ch96 = 1 + [b(:, [2 3]) Y.nb - 1 b(:, 1)]*[1 4 16 48]';

Y = sort_struct(Y, {'ch96'});
[u, ui] = unique(Y.ch96);

%add hairpin covariates
P.covar_dir = cat(2, covar_blocks{:}, 'ref/cod_APOBEC.lawrence.hairpin_covars.align75_filtered/v2', 'ref/cod_c512_covars.align75_filtered/v1');
P.random_effect_cat_covar_index = 'last';

%init the runfile's handle
f = 9;
fs = fopen(sprintf('%s/runs_%d.cmd', basedir, f), 'w');

for x = [ui [ui(2:end) - 1; slength(Y)]]',
  k = x(1); l = x(2);

  %skip if too few mutations for meaningful fit
  %50 is empirically chosen cutoff.
  if nnz(M.mut.count_nb(M.mut.tier == 4 & ismember(M.mut.c512, Y.c512(k:l)), Y.nb(k)) >= 2) < 50, continue; end

  contexts = strjoin(num2cellstr([Y.c512(k:l); c1024(Y.c512(k:l), 2)]), ',');

  %number of patients for this channel
  P.log_exposure = log(full(M.npat(Y.ch96(k))));

  paramspath = sprintf('%s/params/params_%d_%d.mat', basedir, Y.ch96(k), f);
  P.output_name = sprintf('%d-%s_%d', Y.ch96(k), name96{Y.ch96(k)}, f);

  save(paramspath, 'P')

  fprintf(fs, 'mcc/pois_LN_reg_wrapper %s %s c1025 %d %s\n', mafpath, contexts, Y.nb(k), paramspath);
end

% manually add runs for channel 94 (allowing us to assess RXRA)
k = ui(find(u == 94)); l = ui(find(u == 94) + 1) - 1;

contexts = strjoin(num2cellstr([Y.c512(k:l); c1024(Y.c512(k:l), 2)]), ',');

P.log_exposure = log(full(M.npat(Y.ch96(k))));

paramspath = sprintf('%s/params/params_%d_%d.mat', basedir, Y.ch96(k), f);
P.output_name = sprintf('%d-%s_%d', Y.ch96(k), name96{Y.ch96(k)}, 9);

save(paramspath, 'P')

fprintf(fs, 'mcc/pois_LN_reg_wrapper %s %s c1025 %d %s\n', mafpath, contexts, Y.nb(k), paramspath);

% we also need to run this sans hairpin track (since it wasn't run in the first place)

P.covar_dir = cat(2, covar_blocks{:}, 'ref/cod_c512_covars.align75_filtered/v1');

f = 3;

paramspath = sprintf('%s/params/params_%d_%d.mat', basedir, Y.ch96(k), f);
P.output_name = sprintf('%d-%s_%d', Y.ch96(k), name96{Y.ch96(k)}, f);

save(paramspath, 'P')

fprintf(fs, 'mcc/pois_LN_reg_wrapper %s %s c1025 %d %s\n', mafpath, contexts, Y.nb(k), paramspath);

fclose(fs);

% }}}

%}}}

%}}}

%
% hypermutant split subcohorts {{{

clear


X = generate_signature_subcohort_definitions('mutation_data/MC3.align75.ICE_PoN-uniqued.sig_cohorts_NMF/*.M');
X = rmfield(X, 'mafpath');

%update MAFpath to reflect splits
F = [];
F.file = direc('mutation_data/MC3.align75.ICE_PoN-uniqued.mutratesplit_NMF/*');
F = parsein(F, 'file', '.*/(.*)_(hi|lo).M$', {'sig' 'rate'});
F = sort_struct(F, {'sig' 'rate'});

F_hi = reorder_struct(F, strcmp(F.rate, 'hi'));
X.mafpath_hi = mapacross(X.names, F_hi.sig, F_hi.file);
F_lo = reorder_struct(F, strcmp(F.rate, 'lo'));
X.mafpath_lo = mapacross(X.names, F_lo.sig, F_lo.file);

%RC LuT
c1024 = dec2base(0:1023, 4) - 48;
c1024 = [(1:1024)' sum(bsxfun(@times, 3 - c1024(:, [1 3 2 5 4]), [4^4 4^3 4^2 4 1]), 2) + 1];

%base params file
load('LNP_posteriors/params.mat', 'P')

P.covar_dir = {'ref/cod_dnasehs_covars.zerofilled.align75_filtered/v3' ...
               'ref/cod_nuc_covars.align75_filtered/Gm12878' ...
               'ref/cod_mutsig_transformed_covars.align75_filtered/exprmax_transformed' ...
               'ref/cod_mutsig_transformed_covars.align75_filtered/rt_transformed' ...
               'ref/cod_c512_covars.align75_filtered/v1'};
P.random_effect_cat_covar_index = 'last';

P = rmfield(P, 'covar_names');

P.compute_posterior_predictive = 0;
P.compute_marginal_likelihood = 0;
P.init_from_MAP = 0;
P.interval_list = 'ref/target_list.align75_filtered_tier4.txt';
P.tailor_hyperparameters = 1;
P.niter = 3e3;
P.skip_figures = 1;

%generate ch96 names
bc = {'A->C' 'A->G' 'A->T' 'C->A' 'C->G' 'C->T'};
B = {'A' 'C' 'G' 'T'};

tmp = dec2base((1:96) - 1, 4) - 48;
name96 = strcat(B(tmp(:, end - 1) + 1), '(', bc(sum(bsxfun(@times, tmp(:, 1:end - 2), [4 1]), 2) + 1), ')', B(tmp(:, end) + 1))';

%loop over signature subcohorts
for i = 1:slength(X),
  for rates = {'lo' 'hi'},
    rate = rates{1};
    field = ['mafpath_' rate]; 
    mafpath = X.(field){i};
    M = loadM(mafpath);

    basedir = ['LNP_posteriors/signature_subcohorts_mutrate_split/' X.names{i} '_' rate];
    ede(basedir) 
    ede([basedir '/params'])

    P.outdir = [basedir '/output/'];
    ede(P.outdir)

    %open runfile
    f = fopen([basedir '/runs.cmd'], 'w');

    for j = 1:length(X.output{i}),
      Y = [];
      Y.c512 = X.output{i}{j}(:, 1); 
      b = dec2base(Y.c512 - 1, 4, 5) - 48;

      Y.nb = X.output{i}{j}(:, 2);
      Y.ch96 = 1 + [b(:, [2 3]) Y.nb - 1 b(:, 1)]*[1 4 16 48]';

      Y = sort_struct(Y, {'ch96'});
      [~, ui] = unique(Y.ch96);

      for x = [ui [ui(2:end) - 1; slength(Y)]]',
	k = x(1); l = x(2);

	%skip if too few mutations for meaningful fit
	%50 is empirically chosen cutoff.
	if nnz(M.mut.count_nb(M.mut.tier == 4 & ismember(M.mut.c512, Y.c512(k:l)), Y.nb(k)) >= 2) < 50, continue; end

	contexts = strjoin(num2cellstr([Y.c512(k:l); c1024(Y.c512(k:l), 2)]), ',');

	%number of patients for this channel
	P.log_exposure = log(full(M.npat(Y.ch96(k))));

	paramspath = sprintf('%s/params/params_%d.mat', basedir, Y.ch96(k));
	P.output_name = sprintf('%d-%s', Y.ch96(k), name96{Y.ch96(k)});

	save(paramspath, 'P')

	fprintf(f, 'mcc/pois_LN_reg_wrapper %s %s c1025 %d %s\n', mafpath, contexts, Y.nb(k), paramspath);
      end
    end
    fclose(f);
  end
end

%prune output to ensure that only channels with enough mutations for both hi/lo are represented
F = [];
F.params = direc('LNP_posteriors/signature_subcohorts_mutrate_split/*/params/*.mat');
F = parsein(F, 'params', '.*/([A-Z_]+)_(hi|lo).*params_(\d+)\.mat$', {'sig' 'rate' 'ch96'});
F = makeapn(F);
[~, ~, ruj] = unique(F.rate);
[~, ~, suj] = unique(F.sig);
F.uid = mod(F.ch96 + 768*(ruj - 1) + 96*(suj - 1), 768);

hilo_ch_idx = find(accumarray(F.uid, 1) == 2);

F.use = ismember(F.uid, hilo_ch_idx);

save_lines(F.params(F.use), 'LNP_posteriors/signature_subcohorts_mutrate_split/hilo_runs.txt')

%}}}

% }}}

%
% dispatch run files:

% running will take quite some time to finish, so it is recommended that they
% be dispatched to a cluster.

% dispatch the commands output by this bash command:

% find LNP_posteriors/ -name "*.cmd" -exec cat {} +

%}}}

%
% 2. Uniform Poisson {{{
clear

%we can pull covariates from LNP hierarchical runs
F = [];
F.file = direc('LNP_posteriors/pancancer/output/*.mat');
F = parsein(F, 'file', '.*output/(\d+)-(.*)_3.mat', {'ch96' 'context'});
F.ch96 = str2double(F.ch96);
F = sort_struct(F, 'ch96');

ch1536lut = load_struct('ref/1536_LuT.txt');
ch1536lut = makeapn(ch1536lut);
ch1536lut = sparse(ch1536lut.c512, ch1536lut.nbidx, ch1536lut.ch1536);

Beta = NaN(slength(F), 20);

M = cell(slength(F), 1);
pv = cell(slength(F), 1);
XX = cell(slength(F), 1);
n_extra = NaN(slength(F), 1);

pp = parpool(16);
parfor j = 1:slength(F),
  X = load(F.file{j});

  n_extra(j) = nnz(X.M == 0);

  C = [X.C(:, 1:4) full(sparse(1:length(X.C), X.C(:, 5), 1))];
  [BetaP, dBeta, stats] = glmfit(C, X.M, 'poisson', 'constant', 'off');

  %compute p-values  
  pois_p = poisspdf(X.M, exp(C*BetaP));
  pv{j} = pois_p + poisscdf(X.M, exp(C*BetaP), 'upper');
  X.Mu.p_unif = pv{j}(X.M > 0);
  X.Mu.prob_unif = pois_p(X.M > 0);

  %add newbase information
  nb = mod(ceil(F.ch96(j)/16) - 1, 3) + 1;
  X.Mu.ch1536 = full(ch1536lut(X.Mu.c512, nb));

  M{j} = X.Mu;

  X.pois_prob = pois_p;
  X.pois_p = pv{j};
  
  XX{j} = X;
end
pp.delete

idx = ~cellfun(@isempty, XX);

%concat mutation structs and p-value list
L = concat_structs(M(idx));
P = cat(1, pv{:});

%save raw results
for i = find(idx)',
  X = XX{i};
  save(sprintf('pois_reg/output_v1/%d-%s.mat', F.ch96(i), F.context{i}), 'X')
end

%add FDR values to mutation struct
L.q_unif = fdr_jh(L.p_unif, sum(n_extra));

save('pois_reg/output_v1/loci_pvalues.mat', 'L', 'P')

% }}}

%
% Uniform-within-gene {{{
clear

%load territories
load('ref/gene_list.align75_filtered.territories.mat', 'G')

%lookup tables for ch96 -> c32 x 3
tmp = dec2base(0:95, 4) - 48;
c96_32map = [tmp(:, 1:2)*[4 1]' > 2 tmp(:, 3:4)]*[16 4 1]' + 1;
c96_bc = mod(ceil((1:96)/16) - 1, 3)' + 1;

%total number of trimers in exome
N32 = sum(G.terr32);
N96 = N32(c96_32map);

%load mutations
M = loadM('mutation_data/MC3.align75.ICE_PoN.M');

%overall rate of each base change
r96 = histc(M.mut.ch96, 1:96)./N96';
r96_ch = full(sparse(c96_32map, c96_bc, r96));

M.gene.gidx2 = listmap(M.gene.name, G.gene);
M.mut.gidx2 = M.gene.gidx2(M.mut.gene_idx);

%compute Fg's
G.terr96 = G.terr32(:, c96_32map);
G.Fg = accumarray(M.mut.gidx2, 1)./sum(bsxfun(@times, G.terr96, r96'), 2);

save('ref/gene_list.align75_filtered.territories_UWG-Fg.mat', 'G')

%load uniqued mutations
M = loadM('mutation_data/MC3.align75.ICE_PoN-uniqued.M');

M.gene.gidx2 = listmap(M.gene.name, G.gene);
M.mut.gidx2 = M.gene.gidx2(M.mut.gene_idx);

%compute Poisson p-values and individual probabilities
M.mut.p = NaN(slength(M.mut), 3);
M.mut.prob = NaN(slength(M.mut), 3);

%probabilities for each gene having zero mutations each context/change
G.prob0 = NaN(slength(G), 32, 3);

%counts of nonmutated positions in each gene, stratified by context/change
G.n0 = NaN(size(G.prob0));

n_extra = NaN(32, 1);
for i = 1:32,
  idx = M.mut.c32 == i;

  terr = N32(i);
  n_extra(i) = 3*(terr - nnz(idx));

  for j = 1:3,
    ct = M.mut.count_nb(idx, j);
    lams = sum(ct)/terr*G.Fg(M.mut.gidx2(idx));
    M.mut.p(idx, j) = poisscdf(ct, lams, 'upper') + poisspdf(ct, lams);
    M.mut.prob(idx, j) = poisspdf(ct, lams);

    %also save prob for sites in gene with zero mutations
    lams0 = sum(ct)/terr*G.Fg;
    G.prob0(:, i, j) = poisspdf(0, lams0);
    G.n0(:, i, j) = accumarray(M.mut.gidx2(idx), M.mut.count_nb(idx, j) > 0, [slength(G) 1]);
  end
end

G.n0 = bsxfun(@plus, -G.n0, G.terr32);

save('ref/gene_list.align75_filtered.territories_UWG-Fg_zerocounts-and-probs.mat', 'G')

M.mut.q = reshape(fdr_jh(M.mut.p(:), sum(n_extra)), [], 3);

%
%convert to L struct format
C = load_struct('ref/1536_LuT.txt');
C = makeapn(C);
C = sparse(C.c512, C.nbidx, C.ch1536);

Mus = cell(3, 1);
for i = 1:3,
  idx = M.mut.count_nb(:, i) > 0;
  Mus{i} = reorder_struct(M.mut, idx);
  Mus{i}.count = Mus{i}.count_nb(:, i);
  Mus{i}.p = Mus{i}.p(:, i);
  Mus{i}.prob = Mus{i}.prob(:, i);
  Mus{i}.q = Mus{i}.q(:, i);
  Mus{i}.ch1536 = full(C(sub2ind(size(C), Mus{i}.c512, i*ones(slength(Mus{i}), 1))));
end
L = concat_structs(Mus);

save('MC3.align75.ICE_PoN.UWG.results.mat', 'L')

% }}}

%
% 4. Gamma-Poisson (negative binomial) {{{

%
%run wrapper to extract covariates {{{
clear

c1024 = dec2base(0:1023, 4) - 48;
c1024 = [(1:1024)' sum(bsxfun(@times, 3 - c1024(:, [1 3 2 5 4]), [4^4 4^3 4^2 4 1]), 2) + 1];

P = [];
P.annotate_and_exit = 1;
P.pool_basechanges = 0;
P.tier = 1:6;
P.outdir = 'nb/covars';
P.covar_dir = {'ref/cod_dnasehs_covars.zerofilled.align75_filtered/v3' ...
               'ref/cod_nuc_covars.align75_filtered/Gm12878' ...
               'ref/cod_mutsig_transformed_covars.align75_filtered/exprmax_transformed' ...
               'ref/cod_mutsig_transformed_covars.align75_filtered/rt_transformed'};
ede(P.outdir)

pp = parpool(16);
parfor i = 1:512,
  for j = 1:3,
    pois_LN_reg_wrapper('mutation_data/MC3.align75.ICE_PoN-uniqued.M', c1024(i, :), 'c1025', j, P)
  end
end
pp.delete

%}}}

%
%run regressions {{{

clear

%we can pull covariates from LNP hierarchical runs
F = [];
F.file = direc('LNP_posteriors/pancancer/output/*.mat');
F = parsein(F, 'file', '.*output/(\d+)-(.*)_3.mat', {'ch96' 'context'});
F.ch96 = str2double(F.ch96);
F = sort_struct(F, 'ch96');

ch1536lut = load_struct('ref/1536_LuT.txt');
ch1536lut = makeapn(ch1536lut);
ch1536lut = sparse(ch1536lut.c512, ch1536lut.nbidx, ch1536lut.ch1536);

Beta = NaN(slength(F), 20);
alpha = NaN(slength(F), 1);
logL = NaN(slength(F), 1);

M = cell(slength(F), 1);
pv = cell(slength(F), 1);
XX = cell(slength(F), 1);
n_extra = NaN(slength(F), 1);

pp = parpool(16);
parfor j = 1:slength(F),
  X = load(F.file{j});

  n_extra(j) = nnz(X.M == 0);

  if max(X.M) == 1, continue; end

  C = [X.C(:, 1:4) full(sparse(1:length(X.C), X.C(:, 5), 1))];
  r = nbreg(C, X.M);

  Beta(j, :) = r.b;
  alpha(j) = r.alpha;
  logL(j) = r.logL;

  %compute p-values  
  nb_p = nbinpdf(X.M, 1/alpha(j), 1./(1 + alpha(j)*exp(C*Beta(j, :)')));
  pv{j} = nb_p + nbincdf(X.M, 1/alpha(j), 1./(1 + alpha(j)*exp(C*Beta(j, :)')), 'upper');
  X.Mu.p_nb = pv{j}(X.M > 0);
  X.Mu.prob_nb = nb_p(X.M > 0);

  %add newbase information
  nb = mod(ceil(F.ch96(j)/16) - 1, 3) + 1;
  X.Mu.ch1536 = full(ch1536lut(X.Mu.c512, nb));

  M{j} = X.Mu;

  X.nb_prob = nb_p;
  X.nb_p = pv{j};
  
  XX{j} = X;
end
pp.delete

idx = ~cellfun(@isempty, XX);

%concat mutation structs and p-value list
L = concat_structs(M(idx));
P = cat(1, pv{:});

%save raw results
for i = find(idx)',
  X = XX{i};
  save(sprintf('nb/output_v2/%d-%s.mat', F.ch96(i), F.context{i}), 'X')
end

%add FDR values to mutation struct
L.q_nb = fdr_jh(L.p_nb, sum(n_extra));

save('nb/output_v2/loci_pvalues.mat', 'L', 'P')

%}}}

%}}}
