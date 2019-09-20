function effect_bars_v2(eff_obs_known, eff_exp_known, eff_obs_unknown, eff_exp_unknown, M, field, rescue_field, fignum, sil_scale_factor)

if ~exist('sil_scale_factor', 'var'), sil_scale_factor = 1; end

%
% barplot {{{
figure(fignum); clf
set(0, 'CurrentFigure', fignum)
pos = get(fignum, 'Position');
set(fignum, 'Position', [pos(1:2) 310 420])

%known genes
sp1 = subplot('Position', [0.1 0.15 0.38 0.8]);
colormap(sp1, [0.5 0.5 0.5; 0 1 0; 1 0 0])
bar([mean(eff_exp_known); mean(eff_obs_known)], 'stacked', 'EdgeColor', 'none')
line(1*[1 1], quantile(eff_exp_known(:, 1), [0.025 0.975]), 'Color', 'k')
line(1*[1 1], quantile(1 - eff_exp_known(:, 3), [0.025 0.975]), 'Color', 'k')
line(2*[1 1], quantile(eff_obs_known(:, 1), [0.025 0.975]), 'Color', 'k')
line(2*[1 1], quantile(1 - eff_obs_known(:, 3), [0.025 0.975]), 'Color', 'k')
xlim([0 3])
ylim([0 1])

%known false positive overlay
x = mean(eff_obs_known);

scale_factor = sil_scale_factor*eff_obs_known(:, 1)./eff_exp_known(:, 1);
f_mis = eff_exp_known(:, 2).*scale_factor;
f_non = eff_exp_known(:, 3).*scale_factor;

sp3 = axes('Position', [0.58 0.15 0.38 0.8]);
overlay = bar([0 0 0 0 0 0; ...
sil_scale_factor*x(1) (1 - sil_scale_factor)*x(1) mean(f_mis) x(2) - mean(f_mis) mean(f_non) x(3) - mean(f_non)], ...
'stacked', 'EdgeColor', 'none', 'FaceColor', 'k', 'FaceAlpha', 0.2);
sp3.Color = 'none';
sp3.XLim = [0 3];
sp3.YLim = [0 1];
sp3.XTick = [];
overlay(2).FaceColor = 'none';
overlay(4).FaceColor = 'none';
overlay(6).FaceColor = 'none';

line(1*[1 1], x(1) + quantile(f_mis, [0.025 0.975]), 'Color', 'k', 'LineStyle', ':')
line(1*[1 1], x(1) + mean(f_mis) + x(2) - mean(f_mis) + quantile(f_non, [0.025 0.975]), 'Color', 'k', 'LineStyle', ':')

%non-known genes
sp2 = subplot('Position', [0.58 0.15 0.38 0.8]);
colormap(sp2, [0.5 0.5 0.5; 0 1 0; 1 0 0])
bar([mean(eff_exp_unknown); mean(eff_obs_unknown)], 'stacked', 'EdgeColor', 'none')
line(1*[1 1], quantile(eff_exp_unknown(:, 1), [0.025 0.975]), 'Color', 'k')
line(1*[1 1], quantile(1 - eff_exp_unknown(:, 3), [0.025 0.975]), 'Color', 'k')
line(2*[1 1], quantile(eff_obs_unknown(:, 1), [0.025 0.975]), 'Color', 'k')
line(2*[1 1], quantile(1 - eff_obs_unknown(:, 3), [0.025 0.975]), 'Color', 'k')
xlim([0 3])
ylim([0 1])

%non-known false positive overlay
x = mean(eff_obs_unknown);

scale_factor = sil_scale_factor*eff_obs_unknown(:, 1)./eff_exp_unknown(:, 1);
f_mis = eff_exp_unknown(:, 2).*scale_factor;
f_non = eff_exp_unknown(:, 3).*scale_factor;

sp3 = axes('Position', [0.58 0.15 0.38 0.8]);
overlay = bar([0 0 0 0 0 0; ...
sil_scale_factor*x(1) (1 - sil_scale_factor)*x(1) mean(f_mis) x(2) - mean(f_mis) mean(f_non) x(3) - mean(f_non)], ...
'stacked', 'EdgeColor', 'none', 'FaceColor', 'k', 'FaceAlpha', 0.2);
sp3.Color = 'none';
sp3.XLim = [0 3];
sp3.YLim = [0 1];
sp3.XTick = [];
overlay(2).FaceColor = 'none';
overlay(4).FaceColor = 'none';
overlay(6).FaceColor = 'none';

line(2*[1 1], x(1) + quantile(f_mis, [0.025 0.975]), 'Color', 'k', 'LineStyle', ':')
line(2*[1 1], x(1) + mean(f_mis) + x(2) - mean(f_mis) + quantile(f_non, [0.025 0.975]), 'Color', 'k', 'LineStyle', ':')

%labels
sp1.Title.String = 'KCGs';
sp1.Title.FontWeight = 'normal';
sp2.Title.String = 'Non-KCGs';
sp2.Title.FontWeight = 'normal';

sp1.YTickLabels = '';
sp1.YLabel.String = 'Fraction of coding effects';

sp1.XTickLabel = {'Expected' 'Observed'};
sp1.XTickLabelRotation = 45;
sp2.XTickLabel = {'Expected' 'Observed'};
sp2.XTickLabelRotation = 45;

% print some figures we mention in the manuscript text
fprintf('KCG observed/expected sil. frac: %0.2f%%/%0.2f%%\n', 100*mean(eff_obs_known(:, 1)), 100*mean(eff_exp_known(:, 1)));
fprintf('Non-KCG observed/expected sil. frac: %0.2f%%/%0.2f%%\n', 100*mean(eff_obs_unknown(:, 1)), 100*mean(eff_exp_unknown(:, 1)));
fprintf('Fraction of non-KCG nonsilent false positives: %0.2f%%\n', 100*mean((f_mis + f_non)./(1 - eff_obs_unknown(:, 1))));

%}}}

%
% total counts (gene level) barplot {{{

%unique to genes.
M.mut = sort_struct(M.mut, 'gene_idx');
[giu, giui, giuj] = unique(M.mut.gene_idx);

G = [];
G.gene = M.gene.name(giu);
G.tier = M.gene.tier(giu);
G.CV_sig_any = M.mut.CV_sig_any(giui);
G.CV_sig_pancan = M.mut.CV_sig_pancan(giui);
G.sig_dNdS = M.mut.dNdS_sig(giui) | M.mut.dNdS_sig_ttype(giui);
G.oridx = M.mut.oridx(giui);

%minimum q-value for each gene
G.min_q = accumarray(giuj, M.mut.(field), [], @(x) min(x));

if ~isempty(rescue_field),
  demand_field(M.mut, rescue_field)
  G.q_rescue = accumarray(giuj, M.mut.(rescue_field), [], @(x) min(x));
end

%is the most significant hotspot silent?  does it overlap a lncRNA?
G.max_sig_sil = NaN(slength(G), 1);
G.max_sig_sil_lnc = NaN(slength(G), 1);
for x = [giui [giui(2:end) - 1; slength(M.mut)] (1:length(giui))']',
  i = x(1); j = x(2); k = x(3);

  [~, midx] = min(M.mut.(field)(i:j));
  qs = M.mut.(field)(i:j) <= 0.1 & M.mut.effect_idx(i:j) == 5;
  ql = M.mut.(field)(i:j) <= 0.1 & M.mut.effect_idx(i:j) == 5 & M.mut.lncRNAtidx(i:j) > 0;

  G.max_sig_sil(k) = qs(midx) == 1;
  G.max_sig_sil_lnc(k) = ql(midx) == 1;
end

%exclude known false positive mutation calls
G.exclude = accumarray(giuj, M.mut.badidx & M.mut.effect_idx == 5, [], @any) & G.max_sig_sil;

%map in actual dNdS values
X = load('G_dNdS.mat', 'G');
G = mapinto(G, X.G, 'gene', 'name', {'dNdS_hist_fine' 'dNdS_95CI' 'dN_95CI'});

%tumortype specific dNdS values (to avoid filtering genes that are only relevant to a single ttype)
load('G_dNdS_perttype.mat', 'Gtt');
Gtt.good = squeeze((log2(Gtt.dNdS_95CI_ttype(:, 1, :)) >= 0.6 & diff(log2(Gtt.dN_95CI_ttype(:, 1:2, :)), [], 2) < 1.2) | (log2(Gtt.dTdS_95CI_ttype(:, 1, :)) >= 0.4 & diff(log2(Gtt.dT_95CI_ttype(:, 1:2, :)), [], 2) < 1.5));
G.tt_good = ismember(G.gene, Gtt.name(sum(Gtt.good, 2) > 0));

% index genes by (KCG x bad) status
bad_dnds_idx = sum(G.dNdS_hist_fine(:, 401:601), 2)/1e4 > 0.95;
bad_sil_idx = G.max_sig_sil & ~G.max_sig_sil_lnc;
bad_idx = bad_dnds_idx | bad_sil_idx;

known_good = G.tier <= 3 & G.min_q <= 0.1 & ~bad_idx;
known_bad = G.min_q <= 0.1 & bad_idx & G.tier <= 3;

unk_idx = G.tier == 4 & ~G.exclude;
unk_good = unk_idx & G.min_q <= 0.1 & ~bad_idx;
unk_bad = unk_idx & G.min_q <= 0.1 & bad_idx;

% print some figures that we mention in the main text
[~, ~, puj] = unique([M.mut.chr M.mut.pos M.mut.ch1536], 'rows');
fprintf('Total number of sig. hotspots/genes: %d/%d\n', ...
        length(unique(puj(M.mut.(field) <= 0.1))), nnz(G.min_q <= 0.1))
fprintf('Number of KCGs containing hotspots: %d/%d\n', nnz(G.tier <= 3 & G.min_q <= 0.1), nnz(G.tier <= 3))
fprintf('Non-KCG bad fraction: %0.2f%% (%d/%d) (dN/dS: %d, only silent: %d)\n', ...
        100*nnz(unk_bad)/nnz(unk_good | unk_bad), nnz(unk_bad), nnz(unk_good | unk_bad), nnz(bad_dnds_idx & unk_bad), ...
        nnz(bad_sil_idx & unk_bad))

% draw figure
figure(10*fignum); clf
set(0, 'CurrentFigure', 10*fignum)
pause(2)

pos = get(10*fignum, 'Position');
set(10*fignum, 'Position', [pos(1:2) 310 420])

cm = lines;

%
%KCGs {{{

idx = known_good | known_bad;
hist_range_fine = linspace(0, 2, 1e3);
prob_neu = sum(G.dNdS_hist_fine(idx, 401:551), 2)/1e4;
prob_pos = sum(G.dNdS_hist_fine(idx, 701:end), 2)/1e4;
prob_neg = sum(G.dNdS_hist_fine(idx, 1:300), 2)/1e4;

prob_pos(diff(log2(G.dN_95CI(idx, 1:2)), [], 2) > 0.6) = 0;
prob_neu(G.max_sig_sil(idx) & ~G.max_sig_sil_lnc(idx)) = 1;
prob_neu(G.tt_good(idx)) = 0;

[~, si_neu] = sortrows([known_good(idx) prob_neu], [-1 2]);
[~, si_pos] = sortrows([known_good(idx) prob_pos], -[1 2]);

ax_neu1 = axes('Position', [0.1 0.15 0.38 0.8]);
imagesc(ones(size(si_neu)), 'AlphaData', prob_neu(si_neu), 'AlphaDataMapping', 'scaled')
ax_neu1.ALim = [0.5 0.8];
ax_neu1.Color = 'none';
ax_neu1.YDir = 'norm';
colormap(ax_neu1, cm(3, :))

ax_pos1 = axes('Position', [0.1 0.15 0.38 0.8]);
imagesc(ones(size(si_pos)), 'AlphaData', prob_pos(si_pos), 'AlphaDataMapping', 'scaled')
ax_pos1.ALim = [0.68 1];
ax_pos1.Color = 'none';
ax_pos1.YDir = 'norm';
colormap(ax_pos1, cm(1, :))

ax_box1 = axes('Position', [0.1 0.15 0.38 0.8]);
ax_box1.Color = 'none';

%box showing KCGs that are rescued via burden test
if isfield(G, 'q_rescue'),
  n_known_rescued = nnz(G.min_q > 0.1 & G.q_rescue <= 0.1 & G.sig_dNdS & G.tier <= 3);
  n_known_lost = nnz(G.min_q > 0.1 & G.q_rescue <= 0.1 & G.tier <= 3);
  %n_unk_rescued = nnz(G.min_q > 0.1 & G.q_rescue <= 0.1 & G.sig_dNdS & G.tier == 4);

  patch([0.5 1.5 1.5 0.5], [nnz(idx)*[1 1] (n_known_rescued + nnz(idx))*[1 1]], 'k', 'FaceColor', [128 196 229]/255, 'EdgeColor', 'none')
  patch([0.5 1.5 1.5 0.5], [nnz(idx)*[1 1] (n_known_lost + nnz(idx))*[1 1]], 'k', 'FaceColor', 'none', 'LineStyle', ':')
end

patch([0.5 1.5 1.5 0.5], [0 0 nnz(idx)*[1 1]], 'k', 'FaceColor', 'none')

linkaxes([ax_pos1 ax_neu1 ax_box1])
ax_pos1.YLim = [0 1200];
ax_pos1.XLim = [1 - 0.8333 1 + 0.8333];

ax_box1.YTick = 0:200:1200;
ax_box1.YTickLabel = [];

% }}}

%
%non-KCGs {{{

idx = unk_good | unk_bad;
hist_range_fine = linspace(0, 2, 1e3);
prob_neu = sum(G.dNdS_hist_fine(idx, 401:601), 2)/1e4;
prob_pos = sum(G.dNdS_hist_fine(idx, 701:end), 2)/1e4;
prob_neg = sum(G.dNdS_hist_fine(idx, 1:300), 2)/1e4;

prob_pos(diff(log2(G.dN_95CI(idx, 1:2)), [], 2) > 0.6) = 0;
prob_neu(G.max_sig_sil(idx) & ~G.max_sig_sil_lnc(idx)) = 1;

[~, si_neu] = sortrows([unk_good(idx) prob_neu], [-1 2]);
[~, si_pos] = sortrows([unk_good(idx) prob_pos], -[1 2]);

ax_neu2 = axes('Position', [0.58 0.15 0.38 0.8]);
imagesc(ones(size(si_neu)), 'AlphaData', prob_neu(si_neu), 'AlphaDataMapping', 'scaled')
ax_neu2.ALim = [0.5 0.8];
ax_neu2.Color = 'none';
ax_neu2.YDir = 'norm';
colormap(ax_neu2, cm(3, :))

ax_pos2 = axes('Position', [0.58 0.15 0.38 0.8]);
imagesc(ones(size(si_pos)), 'AlphaData', prob_pos(si_pos), 'AlphaDataMapping', 'scaled')
ax_pos2.ALim = [0.68 1];
ax_pos2.Color = 'none';
ax_pos2.YDir = 'norm';
colormap(ax_pos2, cm(1, :))

ax_box2 = axes('Position', [0.58 0.15 0.38 0.8]);
ax_box2.Color = 'none';

%box showing non-KCGs that are rescued via burden test
if isfield(G, 'q_rescue'),
  n_unk_rescued = nnz(G.min_q > 0.1 & G.q_rescue <= 0.1 & G.sig_dNdS & G.tier == 4);
  n_unk_lost = nnz(G.min_q > 0.1 & G.q_rescue <= 0.1 & G.tier == 4);

  patch([0.5 1.5 1.5 0.5], [nnz(idx)*[1 1] (n_unk_rescued + nnz(idx))*[1 1]], 'k', 'FaceColor', [128 196 229]/255, 'EdgeColor', 'none')
  patch([0.5 1.5 1.5 0.5], [nnz(idx)*[1 1] (n_unk_lost + nnz(idx))*[1 1]], 'k', 'FaceColor', 'none', 'LineStyle', ':')
end

patch([0.5 1.5 1.5 0.5], [0 0 nnz(idx)*[1 1]], 'k', 'FaceColor', 'none')

linkaxes([ax_pos2 ax_neu2 ax_box2])
ax_pos2.YLim = [0 1200];
ax_pos2.XLim = [1 - 0.8333 1 + 0.8333];

ax_box2.YTick = 0:200:1200;

% }}}

%global properties
ax_neu1.Box = 'on';
ax_neu2.Box = 'on';

ax_neu1.Title.String = 'KCGs';
ax_neu1.Title.FontWeight = 'normal';
ax_neu2.Title.String = 'Non-KCGs';
ax_neu2.Title.FontWeight = 'normal';

ax_neu2.YAxis.FontSize = 8;

ax_neu1.YTickLabels = '';
ax_neu1.YLabel.String = 'Number of genes with sig. hotspots';

ax_neu1.XTick = [];
ax_pos1.XTick = [];
ax_box1.XTick = [];
ax_neu2.XTick = [];
ax_pos2.XTick = [];
ax_box2.XTick = [];

ax_neu1.YTick = [];
ax_box1.YTick = [];
ax_pos1.YTickLabel = [];
ax_pos2.YTick = [];
ax_box2.YTick = [];

%}}}
