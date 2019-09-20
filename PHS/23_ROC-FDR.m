%load results, map in p-values
load('G_dNdS.mat', 'G')
load('G_dNdS_perttype.mat', 'Gtt');
load('G_dNdS_truepositives.mat', 'G_tp')

load('LNP_posteriors/pancancer/loci_pvalues.mat', 'L');
unifpois = load('pois_reg/output_v1/loci_pvalues.mat', 'L');
uwg = load('MC3.align75.ICE_PoN.UWG.results.mat', 'L');
nb = load('nb/output_v2/loci_pvalues.mat', 'L');

L = multimapinto(L, unifpois.L, {'chr' 'pos' 'ch1536'}, {'chr' 'pos' 'ch1536'}, {'p_unif' 'q_unif'}, {'p_unif' 'q_unif'});
L = multimapinto(L, uwg.L, {'chr' 'pos' 'ch1536'}, {'chr' 'pos' 'ch1536'}, {'p' 'q'}, {'p_uwg' 'q_uwg'});
L = multimapinto(L, nb.L, {'chr' 'pos' 'ch1536'}, {'chr' 'pos' 'ch1536'}, {'p_nb' 'q_nb'}, {'p_nb' 'q_nb'});

Gtmp = [];
[Gtmp.gene, ~, gidx] = unique(L.gene);

%minimum q-value for all sites in each gene
Gtmp.q = accumarray(gidx, L.q(:, 3), [], @min);
Gtmp.q_unif = accumarray(gidx, L.q_unif, [], @min);
Gtmp.q_uwg = accumarray(gidx, L.q_uwg, [], @min);
Gtmp.q_nb = accumarray(gidx, L.q_nb, [], @min);

%number of sites in each gene with q <= 0.1
Gtmp.n_q = accumarray(gidx, L.q(:, 3) <= 0.1);
Gtmp.n_q_unif = accumarray(gidx, L.q_unif <= 0.1);
Gtmp.n_q_uwg = accumarray(gidx, L.q_uwg <= 0.1);
Gtmp.n_q_nb = accumarray(gidx, L.q_nb <= 0.1);

G = mapinto(G, Gtmp, 'name', 'gene', {'q' 'q_unif' 'q_uwg' 'q_nb' 'n_q' 'n_q_unif' 'n_q_uwg' 'n_q_nb'});

%
% load in ttype-specific dN/dS to exclude a few genes from our FP list {{{

Gtt.good = squeeze((log2(Gtt.dNdS_95CI_ttype(:, 1, :)) >= 0.6 & diff(log2(Gtt.dN_95CI_ttype(:, 1:2, :)), [], 2) < 0.9) | (log2(Gtt.dTdS_95CI_ttype(:, 1, :)) >= 0.4 & diff(log2(Gtt.dT_95CI_ttype(:, 1:2, :)), [], 2) < 1.3));

% }}}

%
% map in CGC definitions {{{

load('ref/CGC/v1.1.mat', 'C')

% remove genes not annotated as drivers in tumor types we analyzed (i.e., TCGA tumor types)
C.is_KCG(~C.in_TCGA) = false;

%note that M.mut.tier is not subsequently used, but we update to new KCG definitions for posterity
G_tp.CGC = mapacross(G_tp.name, C.GeneSymbol, double(C.is_KCG));
G_tp.CGC(isnan(G_tp.CGC)) = 0;

G_tp = mapinto(G_tp, Gtmp, 'name', 'gene', {'q' 'q_unif' 'q_uwg' 'q_nb'});

%}}}

%
% Figure 2A, Figure S2B {{{

% define TP set
pos_idx = ismember(G.name, G_tp.name(G_tp.CGC == 1 & G_tp.dMdS_change > log2(1.05) & G_tp.dTdS_95CI(:, 1) < 1.5));

% define FP set -- genes under neutral selection (Prob[0.8 < dN/dS < 1.2] > 0.95)
rng = [find(hist_range_fine > 0.8, 1) find(hist_range_fine > 1.2, 1)];
G.dNdS_p = sum(G.dNdS_hist_fine(:, rng(1):rng(2)), 2)/10000;
neg_idx = G.dNdS_p > 0.95;

%
% alternate TP set based on handpicked genes (Table S3, Figure S2B) {{{

tp_list = {'AKT1' ...
'B2M' ...
'BCOR' ...
'BRAF' ...
'CDKN2A' ...
'CHD4' ...
'CIC' ...
'CTNNB1' ...
'EGFR' 'EP300' ...
'ERBB2' ...
'ERBB3' ...
'FBXW7' ...
'FGFR2' ...
'FGFR3' ...
'GNA11' ...
'GNAQ' ...
'GNAS' ...
'HRAS' ...
'IDH1' ...
'IDH2' ...
'KEAP1' ...
'KIT' ...
'KRAS' ...
'MAP2K1' ...
'MAPK1' ...
'MAX' 'MED12' ...
'MET' ...
'MTOR' ...
'NFE2L2' ...
'NOTCH1' ...
'NRAS' ...
'PIK3CA' 'PIK3CB' ...
'PIK3R1' ...
'POLE' ...
'PPP2R1A' ...
'RAC1' ...
'SF3B1' ...
'SMAD4' ...
'SMARCA4' ...
'SPOP' ...
'TP53' ...
'U2AF1' ...
'ACTB' ...
'CCND1' ...
'CDK4' ...
'ERCC2' ...
'GNA13' ...
'KLF5' ...
'MYCN' ...
'RHOB' ...
'SOX17' ...
'STK19'};

pos_idx_r2 = ismember(G.name, tp_list);

% }}}


X = [];
X.fields = {'q' 'q_nb' 'q_unif' 'q_uwg'}';
X.colors = mm_colormap;
X.roc_legends = gobjects(slength(X), 1);
X.fdr_legends = gobjects(slength(X), 1);
X.roc_legend_labels = {'Log-normal-Poisson' 'Gamma-Poisson' 'Uniform Poisson' 'Uniform-within-gene'}';
q_range = [-2 -1];

ax = plot_roc_curve(G, pos_idx, neg_idx, X, q_range, 10);

% LNP method has FPR of 0 even at q = 0.1
% need to manually add in points on the plot where LNP method's TPR falls

q = G.q;
[~, si] = sort(q);
pos_sort = pos_idx(si);
neg_sort = neg_idx(si);

x = log10(cumsum(neg_sort)/nnz(neg_sort));
y = cumsum(pos_sort)/nnz(pos_sort);

qidx = [find(q(si) >= 0.1, 1) find(q(si) >= 0.01, 1)];

x_dummy = x(find(~isinf(x), 1))*[1 1];
scatter(x_dummy, y(qidx), 'Marker', '*', 'MarkerEdgeColor', X.colors(1, :))
text(x_dummy + 0.05, y(qidx), {'10^{-1}' '10^{-2}'}, 'Color', X.colors(1, :))

print('figures/roc_curve.svg', '-dsvg', '-painters')
print('figures/roc_curve-rev.svg', '-dsvg', '-painters')

% version with handpicked TP set (Figure S2B)
ax = plot_roc_curve(G, pos_idx_r2, neg_idx, X, q_range, 11);

pos_sort = pos_idx_r2(si);

y = cumsum(pos_sort)/nnz(pos_sort);

qidx = [find(q(si) >= 0.1, 1) find(q(si) >= 0.01, 1)];

x_dummy = x(find(~isinf(x), 1))*[1 1];
scatter(x_dummy, y(qidx), 'Marker', '*', 'MarkerEdgeColor', X.colors(1, :))
text(x_dummy + 0.05, y(qidx), {'10^{-1}' '10^{-2}'}, 'Color', X.colors(1, :))

print('figures/roc_curve-alt_TP_set.svg', '-dsvg', '-painters')

% for some reason, Illustator cannot read this directly -- need to import to Safari, export as PDF, then import to Illustrator
% in addition, the normal trick of replacing the fonts with sed fails here, because some of the lines are extremely long.
% do this in Illustrator instead

% figures for manuscript text: how many significant hotspots are in genes belonging to truth sets?

lm = listmap(L.gene, G.name);
L.neg_idx = nansub(neg_idx, lm, false);
L.pos_idx = nansub(pos_idx, lm, false);

% number of hotspots in false positive genes
fprintf('LNP: %d\n', nnz(L.q(:, 3) <= 0.1 & L.neg_idx))
fprintf('NB: %d\n', nnz(L.q_nb <= 0.1 & L.neg_idx))
fprintf('Unif: %d\n', nnz(L.q_unif <= 0.1 & L.neg_idx))
fprintf('UWG: %d\n', nnz(L.q_uwg <= 0.1 & L.neg_idx))

% number of hotspots in true positive genes
fprintf('LNP: %d\n', nnz(L.q(:, 3) <= 0.1 & L.pos_idx))
fprintf('NB: %d\n', nnz(L.q_nb <= 0.1 & L.pos_idx))
fprintf('Unif: %d\n', nnz(L.q_unif <= 0.1 & L.pos_idx))
fprintf('UWG: %d\n', nnz(L.q_uwg <= 0.1 & L.pos_idx))

% }}}

%
% Figure 2B {{{

%load mutations
M = loadM('mutation_data/MC3.align75.ICE_PoN-pvalues_a212692_v2.M');

[~, ui, uj] = unique([M.mut.chr M.mut.pos M.mut.ch1536], 'rows');
Mu = reorder_struct(M.mut, ui);
Mu.count = accumarray(uj, 1);
Mu.gene = M.gene.name(Mu.gene_idx);

lm = listmap(Mu.gene, G.name); nidx = ~isnan(lm); lm = lm(nidx);

Mu.pos_idx = false(slength(Mu), 1);
Mu.pos_idx(nidx) = pos_idx(lm);
Mu.neg_idx = false(slength(Mu), 1);
Mu.neg_idx(nidx) = neg_idx(lm);

Mu.q = Mu.q(:, 3);

figure(50); clf
hold on

for m = 1:slength(X),
  %
  %extract relevant q-values
  q = Mu.(X.fields{m});
  [~, si] = sort(q);
  pos_sort = Mu.pos_idx(si);
  neg_sort = Mu.neg_idx(si);

  %
  %q-value vs. FDR figure
  x = -log10(q(si));
  y = cumsum(neg_sort)./cumsum(ones(size(si)));
  X.fdr_legends(m) = plot(x, y, 'Color', X.colors(m, :));

  y_num = cumsum(neg_sort);
  y_den = cumsum(ones(size(si)));

  x1_idx = find(x < 1, 1);
  y1 = y(x1_idx);
  x2_idx = find(x < 2, 1);
  y2 = y(x2_idx);
  x3_idx = find(x < 3, 1);
  y3 = y(x3_idx);
  fprintf('%0.4f (%d/%d) %0.4f (%d/%d) %0.4f (%d/%d) %s\n', ...
           y1, y_num(x1_idx), y_den(x1_idx), y2, y_num(x2_idx), y_den(x2_idx), y3, y_num(x3_idx), y_den(x3_idx), X.fields{m});

  %add q-value threshold circles
  scatter(x([x1_idx x2_idx x3_idx]), y([x1_idx x2_idx x3_idx]), 'Marker', 'o', 'MarkerEdgeColor', X.colors(m, :))

  %add CI
  cs = cumsum(neg_sort);
  csu = cumsum(ones(size(si))) - cs;

  lns = linspace(0, 8, 200);
  rng = NaN(200, 1);
  for i = 1:200,
    rng(i) = find(lns(i) >= -log10(q(si)), 1);
  end
  rng = unique(flipud(rng));

  beta95ci = NaN(size(rng, 1), 2);
  for x = [rng'; 1:length(rng)],
    i = x(1); j = x(2);

    beta95ci(j, :) = icdf('beta', [0.1 0.9], cs(i), csu(i));
  end

  x = -log10(q(si(rng)));
  idx = ~isinf(x) & ~isnan(x) & ~isnan(beta95ci(:, 1));
  fill([flipud(x(idx)); x(idx)], [flipud(beta95ci(idx, 2)); beta95ci(idx, 1)], X.colors(m, :), 'EdgeColor', 'none', 'FaceAlpha', 0.2)
end

ax = gca;

ax.XLim = [0 8];
ax.YLim = [0 0.08];

ctrl = area(0:0.1:8, 10.^-(0:0.1:8), ax.YLim(2), 'LineStyle', '--', 'FaceColor', 0.9*[1 1 1]);
uistack(ctrl, 'bottom')

ax.XTick = 0:8;
ax.XTickLabel = strsplit(sprintf('10^{-%d} ', ax.XTick), ' ');
ax.XTickLabel{1} = '1';

ax.Box = 'on';
ax.Layer = 'top';

%draw lines at q = 0.1/0.01/0.001
for i = 1:3,
  line(i*[1 1], ylim, 'LineStyle', ':', 'Color', 0.25*[1 1 1])
end

legend([ctrl; X.fdr_legends], 'Empirical FDR Threshold', X.roc_legend_labels{:}, 'Location', 'NorthEast')

% XXX: we will add the second line in postprocessing, since adding it here messes up the
%      vertical dimensions
%xlabel({'FDR under model assumptions' '(q value)'})
xlabel('FDR under model assumptions')

ylabel({'Frac. of hotspots in confident false pos. genes' '(Empirical FDR)'})

print('figures/empirical_FDR_vs_q.svg', '-dsvg')

%fix SVG fonts
% ssed -Ri -e "s/(.*font-family:')SansSerif('.*)/\1Helvetica\2/g" figures/empirical_FDR_vs_q.svg
% ssed -Ri -e "s/(.*font-family:')SansSerif('.*)/\1Helvetica\2/g" figures/empirical_FDR_vs_q-rev.svg

% }}}

%
% top 100 significant non-KCGs by each method -- are we enriched for real hits? {{{

% "nonparametric benchmark" at end of section "Comparative Analyses Confirm Improved Performance of the LNP Model"

hist_range_fine = linspace(0, 2, 1e3);
rng = [find(hist_range_fine > 0.8, 1) find(hist_range_fine > 1.2, 1)];
G.dNdS_p = sum(G.dNdS_hist_fine(:, rng(1):rng(2)), 2)/10000;
G.neg_idx = G.dNdS_p > 0.95;

G.is_KCG = mapacross(G.name, C.GeneSymbol, double(C.is_KCG), 0);

%FP genes by dN/dS
for f = {'q' 'q_nb' 'q_unif' 'q_uwg'},
  H = sort_struct(reorder_struct(G, ~G.is_KCG), f);
  idx = H.dNdS_p(1:100) > 0.95;

  fprintf('%s (%d genes):\n', f{1}, nnz(idx))
  pr(H, {'name' 'q' 'q_nb' 'q_unif' 'q_uwg'}, idx)
  fprintf('\n')
end

%print whole list for each
for f = {'q_nb' 'q_unif' 'q_uwg'},
  H = sort_struct(reorder_struct(G, ~G.is_KCG), f);
  H.n = (1:slength(H))';
  idx = find(H.q > 0.1 & H.n <= 100);

  fprintf('%s (%d genes):\n', f{1}, nnz(idx))
  pr(H, {'n' 'name' 'q' 'q_nb' 'q_unif' 'q_uwg'}, idx)
  fprintf('\n')
end

% }}}
