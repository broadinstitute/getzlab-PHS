M = loadM('mutation_data/MC3.align75.ICE_PoN-pvalues_a212692_v2.M');

%collapse to min. q-value per gene
M.mut = sort_struct(M.mut, 'gene_idx');

[giu, giui, giuj] = unique(M.mut.gene_idx);
G = [];
G.gene = M.gene.name(giu);
G.tier = M.gene.tier(giu);

G.min_q_lnp = accumarray(giuj, M.mut.q(:, 3), [], @(x) min(x));
G.min_q_unif = accumarray(giuj, M.mut.q_unif, [], @(x) min(x));
G.min_q_uwg = accumarray(giuj, M.mut.q_uwg, [], @(x) min(x));
G.min_q_nb = accumarray(giuj, M.mut.q_nb, [], @(x) min(x));

%
% index TP/FP genes based on dN/dS {{{

%load in dN/dS analysis (from 22_gene_dNdS_analysis)
X = load('G_dNdS.mat');

% identify FP genes (neutral dN/dS)
hist_range_fine = linspace(0, 2, 1e3);
rng = [find(hist_range_fine > 0.8, 1) find(hist_range_fine > 1.2, 1)];
X.G.dNdS_p = sum(X.G.dNdS_hist_fine(:, rng(1):rng(2)), 2)/10000;
X.G.bad = X.G.dNdS_p > 0.95 & ~ismember(X.G.name, Gtt.name(sum(Gtt.good, 2) > 0));
G = mapinto(G, X.G, 'gene', 'name', {'bad'}, {'bad_dNdS'});

% identify TP genes (genes whose dM/dS changes by more than 5% when removing hotspots significant
% by UWG model at q < 0.25)
load('G_dNdS_truepositives.mat');
G.true_pos = ismember(G.gene, G_tp.name(G_tp.tier < 4 & G_tp.dMdS_change > log2(1.05) & G_tp.dTdS_95CI(:, 1) < 1.5));

% }}}

%
% plot Figure S2A {{{

r = {'lnp' 'uwg' 'unif'};
c = {'uwg' 'unif' 'nb'};

r = {'nb' 'unif' 'uwg'};
c = {'lnp' 'uwg' 'unif'};

c = {'uwg' 'unif' 'nb'};
r = {'lnp' 'nb' 'unif'};

fullname = [{'uwg' 'unif' 'nb' 'lnp'}; {'Uniform-within-gene' 'Uniform Poisson' 'Gamma-Poisson' 'Log-normal-Poisson'}];

figure(1); clf
set(1, 'Position', [400 0 1380 1140])

ax = gobjects(3, 3);
for i = 1:3, for j = 1:(3 - i + 1),
  ax(i, j) = q_scatter_2(G, c{i}, r{j});

  if i == 1,
    ax(i, j).YLabel.String = mapacross(r{j}, fullname(1, :), fullname(2, :));
  end

  if j == 3 - i + 1,
    ax(i, j).XLabel.String = mapacross(c{i}, fullname(1, :), fullname(2, :));
  end

  ax(i, j).Position = [0.9*(i - 1), 3 - j, 0.85, 0.85]*(1/3) + 0.05*[1 0.85 0 0];
end, end

% version in original draft
print('figures/q_scatter.svg', '-dsvg', '-painters')

% clean SVG
% ssed -Ri -e "s/(.*font-family:')SansSerif('.*)/\1Helvetica\2/g" -e "s/(.*font-family:')mwb_cmsy10('.*)/\1Arial Unicode\2/g" -e "s/>5</>\&#x2264;</g" figures/q_scatter.svg

% }}}
