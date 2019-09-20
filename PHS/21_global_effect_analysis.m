%
% 0. preprocess (load mutations, make O structs) {{{
clear

%
% 0.1. load mutations/p-values  {{{

%load full set of mutations
M = loadM('mutation_data/MC3.align75.ICE_PoN.M');

%RC newbase
M.mut.newbase_idx_rc = M.mut.newbase_idx;
idx = M.mut.c1025 > 512;
M.mut.newbase_idx_rc(idx) = 5 - M.mut.newbase_idx(idx); 

%load p-values
load('LNP_posteriors/pancancer/loci_pvalues.mat', 'L');
gampois = load('nb/output_v2/loci_pvalues.mat', 'L');
unifpois = load('pois_reg/output_v1/loci_pvalues.mat', 'L');
uwg = load('MC3.align75.ICE_PoN.UWG.results.mat', 'L');

%map in p-values
M.mut = multimapinto(M.mut, L, {'chr' 'pos' 'ch1536'}, {'chr' 'pos' 'ch1536'}, {'p' 'q'});
M.mut = multimapinto(M.mut, gampois.L, {'chr' 'pos' 'ch1536'}, {'chr' 'pos' 'ch1536'}, {'p_nb' 'q_nb'});
M.mut = multimapinto(M.mut, unifpois.L, {'chr' 'pos' 'ch1536'}, {'chr' 'pos' 'ch1536'}, {'p_unif' 'q_unif'});
M.mut = multimapinto(M.mut, uwg.L, {'chr' 'pos' 'ch1536'}, {'chr' 'pos' 'ch1536'}, {'p' 'q'}, {'p_uwg' 'q_uwg'});

%
% 0.1.1. index genes recoverable by burden test, events overlapping lncRNAs, or events overlapping TFBSs {{{

M.mut.gene = M.gene.name(M.mut.gene_idx);

%
%index genes whose own dN/dS is confidently high {{{
load('G_dNdS.mat')

G.good = (log2(G.dNdS_95CI(:, 1)) >= 0.25 & diff(log2(G.dN_95CI(:, 1:2)), [], 2) < 0.6) | (log2(G.dTdS_95CI(:, 1)) >= 0.25 & diff(log2(G.dT_95CI(:, 1:2)), [], 2) < 0.75);

M.mut.dNdS_sig = ismember(M.mut.gene, G.name(G.good));

% }}}

%
%index lncRNA overlaps {{{
R = load_refseq('ref/gencode/lncRNA_refSeq.mat');
R.chr = convert_chr(R.chr);
[~, ~, R.uid] = unique(R.id);
cs = cumsum(R.n_exons);
R.range = [[1; cs(1:(end - 1)) + 1] cs];
n_exons = sum(R.n_exons);

%convert to target list.  when mapping, we'll use whichever isoform appears first.
T = [];
T.gene = cell(n_exons, 1);
T.chr = NaN(n_exons, 1);
T.start = NaN(n_exons, 1);
T.end = NaN(n_exons, 1);

for i = 1:slength(R),
  a = R.range(i, 1); b = R.range(i, 2);

  T.gene(a:b) = repmat(R.transcript(i), b - a + 1, 1);
  T.chr(a:b) = R.chr(i);
  T.start(a:b) = max(R.code_start(i), R.exon_starts{i});
  T.end(a:b) = min(R.code_end(i), R.exon_ends{i});
end
T = sort_struct(T, {'chr' 'start' 'end'});

%map mutations to targets
M.mut = sort_struct(M.mut, {'chr' 'pos'});
M.mut.lncRNAtidx = map_mutations_to_targets_fast([M.mut.chr M.mut.pos], [T.chr T.start T.end]);

%}}}

%
% index TFBS overlaps {{{

clear

load('ref/tfbs/wgEncodeRegTfbsClusteredV3.bed.mat', 'T')

T = sort_struct(T, {'chr' 'start' 'end'});

Te = load_struct('ref/gencode/target_list.v1.txt');
Te = makeapn(Te);
Te = sort_struct(Te, {'chr' 'start' 'end'});

[~, cui_TE] = unique(Te.chr);
bdy_TE = [cui_TE [cui_TE(2:end) - 1; slength(Te)]];
[~, cui_T] = unique(T.chr);
bdy_T = [cui_T [cui_T(2:end) - 1; slength(T)]];

% loop over chromosomes; 
TI_ints = cell(24, 1);
pp = parpool(12);
parfor c = 1:24,
  r1 = bdy_TE(c, 1); r2 = bdy_TE(c, 2);
  TE_ints = [Te.start(r1:r2) Te.end(r1:r2)];

  r1 = bdy_T(c, 1); r2 = bdy_T(c, 2);
  T_ints = [T.start(r1:r2) T.end(r1:r2)];

  % filter TF binding sites by score
  n = T.num_exp(r1:r2);
  f_1k = cellfun(@(x) nnz(x > 800), T.exp_scores(r1:r2))./T.num_exp(r1:r2);

  % we keep only sites tested in more than 90 experiments with >70% of scores > 800
  T_ints_cons = T_ints(n > 90 & f_1k > 0.7, :);

  I_ints = interval_list_intersect(TE_ints, T_ints_cons);

  tmp = [];
  tmp.start = I_ints(:, 1);
  tmp.end = I_ints(:, 2);
  tmp.chr = c*ones(slength(tmp), 1);

  TI_ints{c} = tmp;
end

TI_ints = concat_structs(TI_ints);

save('TFBS_overlap_intervals_90-800.mat', 'TI_ints')

% }}}

%}}}

M.mut = rmfield(M.mut, 'gene');
saveM(M, 'mutation_data/MC3.align75.ICE_PoN-pvalues_a212692_v2.M')

%}}}

%
% 0.2. generate O structs {{{

%
% for new definition of KGCs (from CGC March 2018) {{{

%load definition of KCGs
load('ref/CGC/v1.1.mat', 'C')

M.mut.tier = mapacross(M.gene.name(M.mut.gene_idx), C.GeneSymbol, double(C.is_KCG));
M.mut.tier(M.mut.tier ~= 1) = 4;

%
% 0.2.0. for UWG model {{{

%known cancer genes
kidx = M.mut.tier <= 3 & M.mut.q_uwg <= 0.1;

%non-known cancer genes
uidx = M.mut.tier == 4 & M.mut.q_uwg <= 0.1 & ~(M.mut.lncRNAtidx > 0 & M.mut.effect_idx == 5);

O_known = [];
O_known.ef = NaN(512, 3); %observed effects
O_known.nb = NaN(512, 3); %observed newbase
O_known.any = false(512, 1);

O_unknown = O_known;

nbidx = [0 1 2 3; ...
         1 0 2 3; ...
         1 2 0 3; ...
         1 2 3 0];

for j = 1:512,
  chidx = M.mut.c512 == j & ismember(M.mut.effect_idx, [1 3 5]); %& vanidx;
  ch_base = dec2base(j - 1, 4, 5) - 47; ch_base = ch_base(1);

  %if no significant loci in either 
  if ~any((chidx & kidx) | (chidx & uidx)), continue; end

  idx = chidx & kidx;
  if any(idx),
    O_known.ef(j, [2 3 1]) = histc(M.mut.effect_idx(idx), [1 3 5]);
    O_known.nb(j, :) = histc(nbidx(ch_base, M.mut.newbase_idx_rc(idx)), 1:3);
    O_known.any(j) = true;
  end
  idx = chidx & uidx;
  if any(idx),
    O_unknown.ef(j, [2 3 1]) = histc(M.mut.effect_idx(idx), [1 3 5]);
    O_unknown.nb(j, :) = histc(nbidx(ch_base, M.mut.newbase_idx_rc(idx)), 1:3);
    O_unknown.any(j) = true;
  end
end

save('figures/effect_counts-uwg-a212692_CGC.mat', 'O_known', 'O_unknown')

%}}}

%
% 0.2.1. for uniform model {{{

%known cancer genes
kidx = M.mut.tier <= 3 & M.mut.q_unif <= 0.1;

%non-known cancer genes
uidx = M.mut.tier == 4 & M.mut.q_unif <= 0.1 & ~(M.mut.lncRNAtidx > 0 & M.mut.effect_idx == 5);

O_known = [];
O_known.ef = NaN(512, 3); %observed effects
O_known.nb = NaN(512, 3); %observed newbase
O_known.any = false(512, 1);

O_unknown = O_known;

nbidx = [0 1 2 3; ...
         1 0 2 3; ...
         1 2 0 3; ...
         1 2 3 0];

for j = 1:512,
  chidx = M.mut.c512 == j & ismember(M.mut.effect_idx, [1 3 5]); %& vanidx;
  ch_base = dec2base(j - 1, 4, 5) - 47; ch_base = ch_base(1);

  %if no significant loci in either 
  if ~any((chidx & kidx) | (chidx & uidx)), continue; end

  idx = chidx & kidx;
  if any(idx),
    O_known.ef(j, [2 3 1]) = histc(M.mut.effect_idx(idx), [1 3 5]);
    O_known.nb(j, :) = histc(nbidx(ch_base, M.mut.newbase_idx_rc(idx)), 1:3);
    O_known.any(j) = true;
  end
  idx = chidx & uidx;
  if any(idx),
    O_unknown.ef(j, [2 3 1]) = histc(M.mut.effect_idx(idx), [1 3 5]);
    O_unknown.nb(j, :) = histc(nbidx(ch_base, M.mut.newbase_idx_rc(idx)), 1:3);
    O_unknown.any(j) = true;
  end
end

save('figures/effect_counts-poiss-a212692_CGC.mat', 'O_known', 'O_unknown')

%}}}

%
% 0.2.2. for LN-Poiss model {{{

%known cancer genes
kidx = M.mut.tier <= 3 & M.mut.q(:, 3) <= 0.1;

%non-known cancer genes
uidx = M.mut.tier == 4 & M.mut.q(:, 3) <= 0.1 & ~(M.mut.lncRNAtidx > 0 & M.mut.effect_idx == 5);

O_known = [];
O_known.ef = NaN(512, 3); %observed effects
O_known.nb = NaN(512, 3); %observed newbase
O_known.any = false(512, 1);

O_unknown = O_known;

nbidx = [0 1 2 3; ...
         1 0 2 3; ...
         1 2 0 3; ...
         1 2 3 0];

for j = 1:512,
  chidx = M.mut.c512 == j & ismember(M.mut.effect_idx, [1 3 5]); %& vanidx;
  ch_base = dec2base(j - 1, 4, 5) - 47; ch_base = ch_base(1);

  %if no significant loci in either 
  if ~any((chidx & kidx) | (chidx & uidx)), continue; end

  idx = chidx & kidx;
  if any(idx),
    O_known.ef(j, [2 3 1]) = histc(M.mut.effect_idx(idx), [1 3 5]);
    O_known.nb(j, :) = histc(nbidx(ch_base, M.mut.newbase_idx_rc(idx)), 1:3);
    O_known.any(j) = true;
  end
  idx = chidx & uidx;
  if any(idx),
    O_unknown.ef(j, [2 3 1]) = histc(M.mut.effect_idx(idx), [1 3 5]);
    O_unknown.nb(j, :) = histc(nbidx(ch_base, M.mut.newbase_idx_rc(idx)), 1:3);
    O_unknown.any(j) = true;
  end
end

save('figures/effect_counts-LNpoiss-a212692_CGC.mat', 'O_known', 'O_unknown')

%}}}

%
% 0.2.3. for gam-Poiss model {{{

%known cancer genes
kidx = M.mut.tier <= 3 & M.mut.q_nb <= 0.1;

%non-known cancer genes
uidx = M.mut.tier == 4 & M.mut.q_nb <= 0.1 & ~(M.mut.lncRNAtidx > 0 & M.mut.effect_idx == 5);

O_known = [];
O_known.ef = NaN(512, 3); %observed effects
O_known.nb = NaN(512, 3); %observed newbase
O_known.any = false(512, 1);

O_unknown = O_known;

nbidx = [0 1 2 3; ...
         1 0 2 3; ...
         1 2 0 3; ...
         1 2 3 0];

for j = 1:512,
  chidx = M.mut.c512 == j & ismember(M.mut.effect_idx, [1 3 5]); %& vanidx;
  ch_base = dec2base(j - 1, 4, 5) - 47; ch_base = ch_base(1);

  %if no significant loci in either 
  if ~any((chidx & kidx) | (chidx & uidx)), continue; end

  idx = chidx & kidx;
  if any(idx),
    O_known.ef(j, [2 3 1]) = histc(M.mut.effect_idx(idx), [1 3 5]);
    O_known.nb(j, :) = histc(nbidx(ch_base, M.mut.newbase_idx_rc(idx)), 1:3);
    O_known.any(j) = true;
  end
  idx = chidx & uidx;
  if any(idx),
    O_unknown.ef(j, [2 3 1]) = histc(M.mut.effect_idx(idx), [1 3 5]);
    O_unknown.nb(j, :) = histc(nbidx(ch_base, M.mut.newbase_idx_rc(idx)), 1:3);
    O_unknown.any(j) = true;
  end
end

save('figures/effect_counts-gampoiss-a212692_CGC.mat', 'O_known', 'O_unknown')

%}}}

%
% 0.2.4. for all mutations in non-KCGs (control) {{{

%non-known cancer genes
uidx = M.mut.tier == 4;

O_unknown = [];
O_unknown.ef = NaN(512, 3); %observed effects
O_unknown.nb = NaN(512, 3); %observed newbase
O_unknown.any = false(512, 1);

nbidx = [0 1 2 3; ...
         1 0 2 3; ...
         1 2 0 3; ...
         1 2 3 0];

for j = 1:512,
  chidx = M.mut.c512 == j & ismember(M.mut.effect_idx, [1 3 5]); %& vanidx;
  ch_base = dec2base(j - 1, 4, 5) - 47; ch_base = ch_base(1);

  %if no significant loci in either 
  if ~any(chidx & uidx), continue; end

  idx = chidx & uidx;
  if any(idx),
    O_unknown.ef(j, [2 3 1]) = histc(M.mut.effect_idx(idx), [1 3 5]);
    O_unknown.nb(j, :) = histc(nbidx(ch_base, M.mut.newbase_idx_rc(idx)), 1:3);
    O_unknown.any(j) = true;
  end
end

save('figures/effect_counts-control-a212692_CGC.mat', 'O_unknown')

%
% }}}

%}}}

%}}}

%}}}

%
% 1. run DMD MCMC {{{
clear


%load territories
load('ref/pentamer_coding_effect_territories_by_tier.align75_filtered.mat', 'effect_terrs')
terr_known = sum(effect_terrs(:, :, :, 1:3), 4);
terr_unknown = sum(effect_terrs(:, :, :, 4), 4);

%
% run on models {{{

pp = parpool(20);

%Uniform-within-gene
load('figures/effect_counts-uwg-a212692_CGC.mat')
[eff_exp_known eff_obs_known] = run_effect_counts_dir_mult_MCMC(O_known, terr_known);
[eff_exp_unknown eff_obs_unknown] = run_effect_counts_dir_mult_MCMC(O_unknown, terr_unknown);

save('figures/effect_distributions_v2-uwg-a212692_CGC.mat', ...
     'eff_exp_unknown', 'eff_obs_unknown', 'eff_exp_known', 'eff_obs_known')

%Uniform Poisson
load('figures/effect_counts-poiss-a212692_CGC.mat')
[eff_exp_known eff_obs_known] = run_effect_counts_dir_mult_MCMC(O_known, terr_known);
[eff_exp_unknown eff_obs_unknown] = run_effect_counts_dir_mult_MCMC(O_unknown, terr_unknown);

save('figures/effect_distributions_v2-poiss-a212692_CGC.mat', ...
     'eff_exp_unknown', 'eff_obs_unknown', 'eff_exp_known', 'eff_obs_known')

%LN-Poisson
load('figures/effect_counts-LNpoiss-a212692_CGC.mat')
[eff_exp_known eff_obs_known] = run_effect_counts_dir_mult_MCMC(O_known, terr_known);
[eff_exp_unknown eff_obs_unknown] = run_effect_counts_dir_mult_MCMC(O_unknown, terr_unknown);

save('figures/effect_distributions_v2-LNpoiss-a212692_CGC.mat', ...
     'eff_exp_unknown', 'eff_obs_unknown', 'eff_exp_known', 'eff_obs_known')

%gam-Poisson
load('figures/effect_counts-gampoiss-a212692_CGC.mat')
[eff_exp_known eff_obs_known] = run_effect_counts_dir_mult_MCMC(O_known, terr_known);
[eff_exp_unknown eff_obs_unknown] = run_effect_counts_dir_mult_MCMC(O_unknown, terr_unknown);

save('figures/effect_distributions_v2-gampoiss-a212692_CGC.mat', ...
     'eff_exp_unknown', 'eff_obs_unknown', 'eff_exp_known', 'eff_obs_known')

%control (all sites in non-KCGs)
load('figures/effect_counts-control-a212692_CGC.mat')
[eff_exp_unknown eff_obs_unknown] = run_effect_counts_dir_mult_MCMC(O_unknown, terr_unknown);

save('figures/effect_distributions_v2-control-a212692_CGC.mat', 'eff_exp_unknown', 'eff_obs_unknown')

%}}}

%}}}

%
% Figure S1 {{{
load('figures/effect_counts-control-a212692_CGC.mat')

%excise code from run_effect_counts_dir_mult_MCMC.m
O_orig = O_unknown;
terr_orig = terr_unknown;

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

C = load_struct('ref/context32_categs.txt');
C = makeapn(C);
C = reorder_struct(C, unique((dec2base((1:512) - 1, 4) - 48)*[16 1 4 0 0]' + 1, 'stable'));
C.terr_1 = squeeze(terr(:, :, 1));
C.terr_2 = squeeze(terr(:, :, 2));
C.terr_3 = squeeze(terr(:, :, 3));

C = sort_struct(C, 'num');

%make figure
mat = [C.terr_1 C.terr_2 C.terr_3];

%need to get 96 channel frequencies; this is what we will scale the plot by.
M = loadM('mutation_data/MC3.align75.ICE_PoN-pvalues_a212692_v2.M');
n96 = histc(M.mut.ch96, 1:96);

figure(1); clf
b = bar(C.terr_1, 'stacked');

bc = {'A->C' 'A->G' 'A->T' 'C->A' 'C->G' 'C->T'};
B = {'A' 'C' 'G' 'T'};

tmp = dec2base((1:96) - 1, 4) - 48;
name96 = strcat(B(tmp(:, end - 1) + 1), '(', bc(sum(bsxfun(@times, tmp(:, 1:end - 2), [4 1]), 2) + 1), ')', B(tmp(:, end) + 1))';

name96_bch = cellfun(@(x) [x(1) '_' x(end)], name96, 'unif', 0);

cm = [0.5 0.5 0.5; 0 1 0; 1 0 0];

figure(2); clf
ax = gobjects(32, 3);
for j = [1:3:9; 3:3:9],
  rng = j(1):j(2);
  jj = j(2)/3;
  for i = 1:32,
    ax(i, jj) = axes;
    b = bar([mat(i, rng)/sum(mat(i, rng)); 0 0 0], 'stacked', 'BarWidth', 1);
    for f = 1:3, b(f).EdgeColor = 'none'; b(f).FaceColor = cm(f, :); end

    ax(i, jj).XLim = [0.5 1.5];
    ax(i, jj).YLim = [0 1];
  end
end

ax_r = [reshape(ax(1:16, :), [], 1); reshape(ax(17:32, :), [], 1)];

f96 = n96./sum(n96);
f96_cs = cumsum([0; f96]);

for i = 1:96,
  ax_r(i).Position = [f96_cs(i) 0.1 f96(i) 0.8];
  ax_r(i).YAxis.Visible = 'off';

  axes(ax_r(i))
  rectangle('Position', [0.5 0 1 1], 'FaceColor', 'none', 'LineWidth', 0.25)

  %draw label if rectangle is large enough
  ax_r(i).XTick = 1;
  ax_r(i).XAxis.TickLabelInterpreter = 'none';
  ax_r(i).XAxis.FontSize = 8;
  ax_r(i).XTickLabel = name96_bch{i};
  ax_r(i).XTickLabelRotation = 90;
end

lego_colors = [1 1 0;0 0.7 0.7;1 0 0;0.1 0.8 0.1;0 0.2 0.8;0.5 0.3 0.7];
lego_colors = lego_colors([5 4 6 2 3 1], :);

bc_text = {'A>C' 'A>G' 'A>T' 'C>A' 'C>G' 'C>T'};

%draw basechange rectangles
bc_width = sum(reshape(f96, 16, []));
bc_x = f96_cs(1:16:(end - 1));
ax_bc = gobjects(6, 1);
for i = 1:6,
  ax_bc(i) = axes('Position', [bc_x(i) 0.92 bc_width(i) 0.07]);
  ax_bc(i).XAxis.Visible = 'off';
  ax_bc(i).YAxis.Visible = 'off';

  rectangle('Position', [0 0 1 1], 'FaceColor', lego_colors(i, :), 'LineWidth', 0.25)

  text(0.5, 0.5, bc_text{i}, 'HorizontalAlignment', 'center', 'Interpreter', 'none', 'FontWeight', 'bold')
end

print('figures/coding_effect_breakdown_by_channel.svg', '-dsvg')

%postprocess SVG
% ssed -Ri -e "s/(.*font-family:')SansSerif('.*)/\1Helvetica\2/g" figures/coding_effect_breakdown_by_channel.svg

%channels incapable of generating various coding effects (in figure legend):

%synonymous
% A->C in A_T
% A->T in A_T
% C->A in A_T
% C->G in A_T

%nonsense
for i = num2cellstr(1:3)',
  idx = find(C.(['terr_' i{1}])(:, 3) == 0);
  C.name(idx)
end

% A in A_A -> C
% A in A_C -> C
% A in A_G -> C
% A in A_T -> C
% 
% C in C_A -> G
% C in C_C -> G
% C in C_G -> G
% C in G_A -> G
% C in G_C -> G
% C in G_G -> G
% C in T_C -> G
% C in T_G -> G
% 
% C in A_T -> T
% C in C_T -> T
% C in G_T -> T
% C in T_T -> T

%}}}

%
% Figure 1 {{{
clear


% load hotspots/p-values
M = loadM('mutation_data/MC3.align75.ICE_PoN-pvalues_a212692_v2.M');

%get CGC definitions
load('ref/CGC/v1.1.mat', 'C')

% remove genes not annotated as drivers in tumor types we analyzed (i.e., TCGA tumor types)
C.is_KCG(~C.in_TCGA) = false;

%note that M.mut.tier is not subsequently used, but we update to new KCG definitions for posterity
M.mut.tier = mapacross(M.gene.name(M.mut.gene_idx), C.GeneSymbol, double(C.is_KCG));
M.mut.tier(M.mut.tier ~= 1) = 4;

%M.gene is used to infer tier in effect_bars_v2
M.gene.tier = mapacross(M.gene.name, C.GeneSymbol, double(C.is_KCG));
M.gene.tier(M.gene.tier ~= 1) = 4;
M.gene.hallmark = mapacross(M.gene.name, C.GeneSymbol, double(strcmp(C.Hallmark, 'Yes')));

% draw figures

load('figures/effect_distributions_v2-uwg-a212692_CGC.mat')
effect_bars_v3(eff_obs_known, eff_exp_known, eff_obs_unknown, eff_exp_unknown, M, 'q_uwg', [], 1, 0.9)
figure(1)
print('figures/effect_distributions_v3-uwg-a212692_CGC_90-sil-neu.eps', '-depsc', '-painters')
figure(10)
print('figures/effect_counts_v3-uwg-a212692_CGC_90-sil-neu.svg', '-dsvg')

load('figures/effect_distributions_v2-poiss-a212692_CGC.mat')
effect_bars_v3(eff_obs_known, eff_exp_known, eff_obs_unknown, eff_exp_unknown, M, 'q_unif', [], 2, 0.9)
figure(2)
print('figures/effect_distributions_v3-unif-a212692_CGC_90-sil-neu.eps', '-depsc', '-painters')
figure(20)
print('figures/effect_counts_v3-unif-a212692_CGC_90-sil-neu.svg', '-dsvg')

load('figures/effect_distributions_v2-gampoiss-a212692_CGC.mat')
effect_bars_v3(eff_obs_known, eff_exp_known, eff_obs_unknown, eff_exp_unknown, M, 'q_nb', 'q_uwg', 3, 0.9)
figure(3)
print('figures/effect_distributions_v3-gampois-a212692_CGC_90-sil-neu.eps', '-depsc', '-painters')
figure(30)
print('figures/effect_counts_v3-gampois-a212692_CGC_90-sil-neu.svg', '-dsvg')

load('figures/effect_distributions_v2-LNpoiss-a212692_CGC.mat')
M.mut.q_lnp_mean = M.mut.q(:, 3);
effect_bars_v3(eff_obs_known, eff_exp_known, eff_obs_unknown, eff_exp_unknown, M, 'q_lnp_mean', 'q_uwg', 4, 0.9)
figure(4)
print('figures/effect_distributions_v3-lnpois-a212692_CGC_90-sil-neu.eps', '-depsc', '-painters')
figure(40)
print('figures/effect_counts_v3-lnpois-a212692_CGC_90-sil-neu.svg', '-dsvg')

%postprocess SVGs

find figures -maxdepth 1 -name "effect_counts_v3-*-a212692_CGC_90-sil-neu.svg" | while read -r i; do
ssed -Ri -e "s/(.*font-family:')SansSerif('.*)/\1Helvetica\2/g" $i
done

%}}}
