%
% Figure 3A {{{

%
% 1. local sampling {{{

%
% 1.1. LNP {{{
clear

F = [];
F.file = direc('LNP_posteriors/pancancer/output/*.mat');
F = parsein(F, 'file', '.*output/(\d+)-(.*)_3.mat', {'ch96' 'context'});
F.ch96 = str2double(F.ch96);
F = sort_struct(F, 'ch96');

full_samps = cell(slength(F), 1);

pp = parpool(20);
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

save('figures/samps_LNP.mat', 'full_samps', '-v7.3')

%}}}

%
% 1.2. Unif. Poiss. {{{
clear

F = [];
F.file = direc('pois_reg/output/*.mat');
F = parsein(F, 'file', '.*output/(\d+)-(.*).mat', {'ch96' 'context'});
F.ch96 = str2double(F.ch96);
F = reorder_struct(F, ~isnan(F.ch96));
F = sort_struct(F, 'ch96');

full_samps = cell(slength(F), 1);

pp = parpool(20);
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

save('figures/samps_UP.mat', 'full_samps', '-v7.3')

%}}}

%
% 1.3. UWG {{{
clear

load('ref/gene_list.align75_filtered.territories_UWG-Fg.mat', 'G')
N32 = sum(G.terr32);

M = loadM('mutation_data/MC3.align75.ICE_PoN-uniqued.M');
M.gene.gidx2 = listmap(M.gene.name, G.gene);
M.mut.gidx2 = M.gene.gidx2(M.mut.gene_idx);

%get full set of lambdas and their counts
lams = NaN(32, 3, slength(G));
for i = 1:32,
  idx = M.mut.c32 == i;

  terr = N32(i);

  for j = 1:3,
    ct = M.mut.count_nb(idx, j);
    lams(i, j, :) = sum(ct)/terr*G.Fg;
  end
end

terr_rep = permute(repmat(G.terr32, 1, 1, 3), [2 3 1]);

u = unique(G.terr32)';

%generate samples
full_samps = cell(100, 1);

pp = parpool(20);
parfor k = 1:100, 
  samps = cell(length(u), 1);
  for x = [u; 1:length(u)]
    i = x(1); j = x(2); 

    idx = terr_rep == i;

    samps{j} = poissrnd(repmat(lams(idx), i, 1));
  end
  full_samps{k} = cat(1, samps{:});
end
pp.delete

%concatenate and save
full_samps = sparse(cat(2, full_samps{:}));

save('figures/samps_UWG.mat', 'full_samps', '-v7.3')

%}}}

%
% 1.4. Gamma-Poiss. {{{
clear

F = [];
F.file = direc('LNP_posteriors/pancancer/output/*.mat');
F = parsein(F, 'file', '.*output/(\d+)-(.*)_3.mat', {'ch96' 'context'});
F.ch96 = str2double(F.ch96);
F = sort_struct(F, 'ch96');

full_samps = cell(slength(F), 1);

pp = parpool(20);
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

save('figures/samps_GP.mat', 'full_samps', '-v7.3')

%}}}

%}}}

%
% 1.1. distributed sampling {{{

%to round out the bottoms of the histograms (especially when computing fractions), we need more
%samples.  not feasible to run this locally; need to run distributed.

%code of each the previous blocks broken out into their own files, which are then
%compiled and dispatched to a cluster

%
% 1.1.1. compile {{{

% cd misc_functions/model_sims/
% find . -name "*.m" | while read -r i; do
% mcc -m -R '-singleCompThread' -d mcc \
% -I ../../funcs \
% $i &
% done

%}}}

%
% 1.1.2. dispatch in parallel {{{

% dispatch the commands output by this bash loop:

% for i in run_{UWG,GP,LNP,UP}; do
% 	for j in {1..10}; do
% 		echo misc_functions/model_sims/mcc/$i $j
% 	done
% done

% }}}

%}}}

%
% 2. draw figure {{{
clear

%
% 2.0. load exome territories (to get observed counts) {{{
%      we'll also index silent positions in simulated vector

F = [];
F.file = direc('LNP_posteriors/pancancer/output/*.mat');
F = parsein(F, 'file', '.*output/(\d+)-(.*)_3.mat', {'ch96' 'context'});
F.ch96 = str2double(F.ch96);
F = sort_struct(F, 'ch96');

%
%index silent mutations/territory {{{

%load covariate matrix to index territories
cmats = cell(24, 1);
for c = 1:24,
  load(sprintf('ref/cod_c512_covars.align75_filtered/v1/chr%d.mat', c), 'covar_mat');
  cmats{c} = covar_mat; 
end

%add ch96 information to 1536 LuT
%this allows mapping ch96 -> c1024 (to pick out relevant sites in the covariate file used to index territory)
lut1536 = load_struct('ref/1536_LuT.txt');
lut1536 = makeapn(lut1536 );
lut1536.c32 = (dec2base(lut1536.c512 - 1, 4) - 48)*[16 1 4 0 0]' + 1;

tmp = dec2base(0:95, 4) - 48;
lut96 = [];
lut96.c96_32map = [tmp(:, 1:2)*[4 1]' > 2 tmp(:, 3:4)]*[16 4 1]' + 1;
lut96.c96_bc = mod(ceil((1:96)/16) - 1, 3)' + 1;
lut96.ch96 = (1:96)';

lut1536 = multimapinto(lut1536, lut96, {'c32' 'nbidx'}, {'c96_32map' 'c96_bc'}, {'ch96'});

save_struct(lut1536, 'ref/1536_LuT_v2.txt')

c1024 = dec2base(0:1023, 4) - 48;
c1024 = [(1:1024)' sum(bsxfun(@times, 3 - c1024(:, [1 3 2 5 4]), [4^4 4^3 4^2 4 1]), 2) + 1];

%load full set of mutations to allow for easy effect lookup
M = loadM('mutation_data/MC3.align75.ICE_PoN-pvalues_a212692.M');

%reverse effect for bases that need RCing
efx = dec2base(0:26, 3) - 47;
efxrc = efx(:, [3 2 1]);
efxrcidx = [1; sum(bsxfun(@times, efxrc - 1, [9 3 1]), 2) + 2; 29];

%for each channel, mask synonymous mutations/territory
syn_mask = cell(96, 1);

%loop over each channel's results
clear X

pp = parpool(20);
parfor j = 1:96,
  %open C+E FWB to index potential effects of nonmutated positions
  fwb = org.broadinstitute.cga.tools.seq.FixedWidthBinary('ref/gencode/c65e29/all.fwb');

  X = load(F.file{j});

  X.Mu.ch96 = j*ones(slength(X.Mu), 1);
  X.Mu.effect_idx = M.mut.effect_idx(multimap(X.Mu, M.mut, {'chr' 'pos' 'ch96'}));

  c1024idx = nonzeros(c1024(lut1536.c512(lut1536.ch96 == j), :));
  nbidx = lut1536.nbidx(find(lut1536.ch96 == j, 1));
  syn_idx = find(efx(:, nbidx) == 1) + 1;

  %mask of synonymous mutations/territory for this channel
  syn_mask{j} = false(size(X.M));

  %subtract mutations from territory; look up which positions in territory are capable of generating synonymous events
  %need to go chromosome-by-chromosome because of how covariates are organized
  [cu, cui] = unique(X.Mu.chr);

  %XXX: this code assumes all chromosomes contain mutations of all 96 channels.
  %     this is an OK assumption for our dataset, but may not work for others.

  offset = 0;

  for x = [cui [cui(2:end) - 1; slength(X.Mu)] cu]',
    i = x(1); k = x(2); c = x(3);

    tmp = cmats{c}(:, c1024idx);
    cm = tmp(:, 1);
    for q = 2:size(tmp, 2),
      cm = cm + tmp(:, q);
    end

    zidx = cm > 0;
    zidx(X.Mu.pos(i:k)) = 0;

    [a, b] = find(tmp(zidx, :));

    [~, si] = sort(a);
    b = b(si);

    con = c1024idx(b);
    eff = mod(double(fwb.get(c, find(zidx))), 29) + 1;

    %flip effect at sites that would also be RC'd
    idx = con > 512; 
    eff(idx) = efxrcidx(eff(idx));

    %finally, mask all sites that would yield a silent mutation given this basechange
    syn_mut_idx = X.Mu.effect_idx(i:k) == 5;

    syn_mask{j}(find(syn_mut_idx) + offset) = true;
    syn_mask{j}(find(ismember(eff, syn_idx)) + k - i + 1 + offset) = true;

    offset = offset + nnz(cm > 0);
  end

  fwb.close();
end
pp.delete

syn_mask = cat(1, syn_mask{:});

save('figures/syn_mask.mat', 'syn_mask')

%}}}

%
% get observed counts {{{

count_hist = NaN(96, 1000);
count_sil_hist = NaN(96, 1000);

terr = NaN(96, 16);
terr_sil = NaN(96, 16);

offset = 0;

for j = 1:slength(F),
  X = load(F.file{j});

  syn_mask_ch = syn_mask((1:size(X.M, 1)) + offset);

  count_hist(j, :) = accumarray(1 + X.M, 1, [1000 1]);
  count_sil_hist(j, :) = accumarray(1 + X.M(syn_mask_ch), 1, [1000 1]);

  terr(j, :) = histc(X.C(:, 5), 1:16);
  terr_sil(j, :) = histc(X.C(syn_mask_ch, 5), 1:16);

  offset = offset + size(X.M, 1);

  fprintf('%d ', j)
end
tot_terr = fullsum(terr);
tot_terr_sil = fullsum(terr_sil);

terr_rng = cumsum(sum(terr, 2)); terr_rng = [[1; terr_rng(1:(end - 1)) + 1] terr_rng];

obs_counts = sum(count_hist);
obs_counts_sil = sum(count_sil_hist);

%}}}

%}}}

%
% 2.1. load in expected counts from each model {{{

%note that for each method, we have 11 total simulated runs: 
% 1 initial run (in the figures directory) and 10 supplementary runs (in figures/model_sims/)

fn_tmp = cell(10, 1);
for i = 1:10,
  fn_tmp{i} = strcat('figures/model_sims/samps_', {'UWG'; 'UP'; 'GP'; 'LNP'}, sprintf('_%d.mat', i));
end


cols = mm_colormap;

X = [];
X.filenames = [strcat('figures/samps_', {'UWG'; 'UP'; 'GP'; 'LNP'}, '.mat') fn_tmp{:}];
X.cols = cols([4 3 2 1], :);
X.h = cell(slength(X), 1);
X.h_syn = cell(slength(X), 1);
X.rng = cell(slength(X), 1);

pp = parpool(10);
for m = 1:slength(X),
  %initialize containers for storing histograms
  h = cell(size(X.filenames, 2), 1);
  h_syn = cell(size(X.filenames, 2), 1);

  %we assume no simulation generated more than 999 mutations at any position 
  rng = 0:999;

  %loop over runs
  for f = 1:size(X.filenames, 2),
    load(X.filenames{m, f}, 'full_samps');

    full_samps_syn = full_samps(syn_mask, :);

    h_tmp = NaN(size(full_samps, 2), length(rng));
    h_syn_tmp = NaN(size(full_samps, 2), length(rng));

    parfor i = 1:size(full_samps, 2),
      h_tmp(i, :) = histc(full(full_samps(:, i)), rng);
      h_syn_tmp(i, :) = histc(full(full_samps_syn(:, i)), rng);
    end

    h{f} = h_tmp;
    h_syn{f} = h_syn_tmp;

    fprintf('%d,%d ', f, m)
  end
  fprintf('\n')

  %concatenate the 11 runs
  X.h{m} = cat(1, h{:});
  X.h_syn{m} = cat(1, h_syn{:});
end
pp.delete

save('figures/model_sims_histograms.mat', 'X')

%}}}

%
% 2.2. plot Figure 3A {{{


[sp1 sp2] = histogram_plot(X.h_syn, obs_counts_sil, nnz(syn_mask), X.cols, 2);
sp1.YLabel.String = 'Fraction of silent exonic positions';
sp1.XLim = [-0.5 15.5];
sp1.YLim = [-10 0];
sp2.YLim = [-0.5 3];
sp1.YTickLabel = strsplit(sprintf('10^{%d} ', sp1.YTick), ' ');
sp1.YTickLabel{end - 1} = sprintf('%0.2f', 10^max(log10(obs_counts_sil/sum(obs_counts_sil))));

print('figures/fig3a.eps', '-depsc', '-painters')

%}}}

%}}}

%}}}

%
% Figure 3B {{{
clear

%
% 1. make base R struct from LNP results {{{
X = load('LNP_posteriors/pancancer/output/3-A(A->C)G_3.mat');

R = [];
R.mut = [];
R.mut.count = X.M;
R.mut.covars = X.C;
f = NaN(size(X.M));
R.mut.chr = f;
R.mut.pos = f;
R.mut.gene_idx = f;

R.mut.tmp = f;
idx = R.mut.count > 0;
R.mut.tmp(idx) = 1:nnz(idx);

X.Mu.tmp = (1:slength(X.Mu))';

R.mut = mapinto(R.mut, X.Mu, 'tmp', {'chr' 'pos' 'gene_idx'});
R.mut = rmfield(R.mut, 'tmp');

R.gene = load2('mutation_data/MC3.M/gene');

%}}}

%
% 2. map in predicted fractions, q-values, and dN/dS status for each method {{{

M = loadM('mutation_data/MC3.align75.ICE_PoN-pvalues_a212692_v2.M');
Mu_ch = reorder_struct(M.mut, M.mut.ch96 == 3);
Mu_ch.q_LNP = Mu_ch.q(:, 3);

R.mut = multimapinto(R.mut, Mu_ch, {'chr' 'pos'}, {'chr' 'pos'}, {'q_LNP' 'q_nb' 'q_unif' 'q_uwg'}, {'q_LNP' 'q_NB' 'q_UP' 'q_UWG'});

% 1. LNP
R.mut.p_LNP = X.post_pred(:, 3);

% 2. NB
load('nb/output_v2/3-A(A->C)G.mat', 'X');
R.mut.p_NB = X.nb_prob;

% 3. UP
load('pois_reg/output_v1/3-A(A->C)G.mat', 'X');
R.mut.p_UP = X.pois_prob;

% 4. UWG

%for some reason the struct only contains ch1536 information; map to ch96
lut1536 = load_struct('ref/1536_LuT_v2.txt');
lut1536 = makeapn(lut1536);

load('MC3.align75.ICE_PoN.UWG.results.mat', 'L');
L = reorder_struct(L, ismember(L.ch1536, lut1536.ch1536(lut1536.ch96 == 3)));

R.mut = multimapinto(R.mut, L, {'chr' 'pos'}, {'chr' 'pos'}, {'prob'}, {'p_UWG'});

% }}}

%
% 3. plot Figure 3B {{{


F = [];
F.field = {'LNP' 'NB' 'UP' 'UWG'}';
F.cols = mm_colormap;
F.cols_num = [0 0 1; 0 0 0; 0 1 0; 1 0 0];
F.cols_light = [0.4 0.4 1; 0.6 0.6 0.6; 0.4 1 0.4; 1 0.4 0.4];

%observed fractions (weighted geometric mean across pentamers)
h = full(sparse(R.mut.count + 1, R.mut.covars(:, 5), 1));
h = bsxfun(@rdivide, h, sum(h));
hc = full(sparse(R.mut.count + 1, R.mut.covars(:, 5), 1));

wgm = NaN(size(h, 1), 1);
for i = 1:size(h, 1),
  idx = h(i, :) > 0;
  wgm(i) = sum(bsxfun(@times, log10(h(i, idx)), hc(i, idx)))/sum(hc(i, idx), 2);
end

%add jitter to x-axis
rng(2345);
x_offset = 0.5*rand(size(R.mut.count, 1), 4) - 0.25;

figure(3); clf

%define subplots
sp = NaN(16, 1);
ax = gobjects(16, 1);
for i = 1:15,
  sp(i) = subplot('Position', [0.0588*(i + 1) 0.12 0.05 0.85]);
  ax(i) = gca;
  hold(ax(i), 'on')
end

%populate subplots
for j = 1:4,
  x = x_offset(:, j) + j - 0.5;
  y = log10(R.mut.(['p_' F.field{j}]));

  for i = 1:15,
    subplot(sp(i))

    cidx = R.mut.count == (i - 1);
    nidx = R.mut.q_NB < 0.1 & R.mut.count >= 7;

    scatter(x(cidx), y(cidx), 'Marker', 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', F.cols(j, :), 'SizeData', 15)
    scatter(x(cidx & nidx), y(cidx & nidx), 'Marker', 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'none', 'SizeData', 25)
  end
end

%axes postprocessing
for i = 1:15,
  ax(i).XLim = [0 4];
  ax(i).YLim = [-40 0];

  ax(i).XTick = [];
  ax(i).Box = 'on';

  ax(i).XLabel.String = i - 1;

  ax(i).XAxis.MinorTickValues = 1:3;
  ax(i).XMinorGrid = 'on';

  if i > 1,
    ax(i).YTickLabel = {};
  else
    ax(i).YTick = -40:5:0;
    ax(i).YTickLabel = strsplit(sprintf('10^{%d} ', ax(i).YTick), ' ');
  end
end

ax(1).YLabel.String = 'Prob(# mutations at each site)';

%draw sig. q patches/observed fractions
q_cut = 4.88e-7;

%draw observed fractions/q cutoff patches
for i = 1:15,
  subplot(sp(i))

  %observed fractions
  if ~isnan(wgm(i)),
    line(xlim, wgm(i)*[1 1], 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', ':')
  end

  %patch for q > 0.1
  p = patch([xlim fliplr(xlim)], [0 0 [1 1]*log10(q_cut)], [0.5 0.5 0.5], 'EdgeColor', 'none', 'FaceAlpha', 0.25);
  ax(i).Children = ax(i).Children([2:end 1]);
end

%draw lines between loci we want to highlight
nidx = R.mut.q_NB < 0.1 & R.mut.count >= 7;
[~, ~, uj] = unique(R.mut.gene_idx(nidx));
x = bsxfun(@plus, x_offset(nidx, :), 1:4) - 0.5;
y = log10(struct2matrix(reorder_struct(R.mut, nidx), {'p_LNP' 'p_NB' 'p_UP' 'p_UWG'}));

gidx = nidx & R.mut.q_LNP < 0.1; gidx = gidx(nidx) + 1;

nnidx = find(nidx);

cols = [1 0 0; 0 0 1];

for i = 1:size(x, 1),
  plot_n = R.mut.count(nnidx(i)) + 1;
  subplot(sp(plot_n))
  
  [y_s, si] = sort(y(i, :));
  line(x(i, si), y_s, 'Color', [cols(gidx(i), :) 0.5], 'LineWidth', 1.5)

  ax(plot_n).Children = ax(plot_n).Children([2:(end - 1) 1 end]);
end

print('figures/fig3b.eps', '-painters', '-depsc')

%}}}

%}}}
