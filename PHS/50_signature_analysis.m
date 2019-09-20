%
% Figure 4A {{{

clear

watercolor_plot(1, 'LNP_posteriors/signature_subcohorts/*/output/*.mat')

% this figure cannot render as a decent vector graphic due to a bug in MATLAB's graphics engine
% thus, we must save as a very high-res raster image.
print('figures/watercolor.png', '-dpng', '-r600')

%}}}

%
% Figure 4B {{{

% note that the most significant bit is flipped, i.e. it correspond to *no* 5-mer categorical
% random effects.

%   !5m fin crs
% 0           
% 1         x
% 2     x
% 3     x   x
% 4 x         
% 5 x       x
% 6 x   x
% 7 x   x   x

%   5m fin crs
% 0 x         
% 1 x       x
% 2 x   x
% 3 x   x   x
% 4           
% 5         x
% 6     x
% 7     x   x

% we thus invert this to make ordering more intuitive.

clear

F = [];
F.file = direc('LNP_posteriors/signature_subcohorts/*/output/*.mat');
F = parsein(F, 'file', '.*/(.*?)/output/(\d+)-([ACGT]\(.->.\)[ACGT])_(\d).*$', {'sig' 'ch96' 'context' 'covars'});
F = makeapn(F);
F = reorder_struct(F, ~(strcmp(F.sig, 'PANCAN') | F.covars > 7 | strcmp(F.sig, ''))); %exclude XRseq UV (covars 8,9)

%invert most significant bit
tmp = dec2base(0:7, 2) - 48; tmp(:, 1) = xor(tmp(:, 1), 1);
F.covars = mapacross(F.covars, 0:7, tmp*[4 2 1]');

F = sort_struct(F, {'sig' 'ch96' 'covars'});

[~, sui] = unique(F.sig);
V = cell(length(sui), 1);

%define binary coefficient matrix
categs = cell(8, 1);
for i = 0:7,
  categs{i + 1} = repmat(dec2base(i, 2, 3) - 48, 2501, 1);
end
categs = repmat(cat(1, categs{:}), 16, 1);

%
% 1. load in data {{{

%for each signature
for x = [sui [sui(2:end) - 1; slength(F)] (1:length(sui))']', 
  Fsig = reorder_struct(F, x(1):x(2));

  %within each trimer
  U = [];
  [U.ch96, cui] = unique(Fsig.ch96); 
  U.sig = repmat(F.sig(x(1)), slength(U), 1);
  U.coeffs_sig2_tri = NaN(slength(U), 4);

  %names of pentamers within this trimer (convert uniqued categories into pentamers)
  U.c512 = NaN(slength(U), 16);

  %sigma_tot^2 == sigma_0^2 + sigma_j^2 (for no covariates/all covariates)
  U.sig2_tot_actual = NaN(slength(U), 16, 2);

  %histogram sigma_tot^2 with all covariates
  U.sig2_tot_hist = NaN(slength(U), 16, 101);

  U.sig2_tot_CI95 = NaN(slength(U), 16, 2);

  for y = [cui [cui(2:end) - 1; slength(Fsig)] (1:length(cui))']', 
    s2s = NaN(2501, 8, 16);
    for z = [y(1):y(2); 1:8],
      i = z(1); j = z(2);

      X = load(Fsig.file{i});

      s2_tots = bsxfun(@plus, 1./X.tau(1:(end - 1), 500:end), 1./X.tau(end, 500:end))';
      if Fsig.covars(i) == 0, 
	%get names of pentamers within this trimer
	%because we didn't do any downsampling for the trimer runs, we know that all 16 will
	%be present, in order.
	U.c512(y(3), :) = unique(X.Mu.c512);

	U.sig2_tot_actual(y(3), :, 1) = mean(s2_tots);
      end
      if Fsig.covars(i) == 7,
	U.sig2_tot_actual(y(3), :, 2) = mean(s2_tots);
	h = histc(s2_tots, 0:0.05:5)'; h = h/2501;
	U.sig2_tot_hist(y(3), :, :) = h;
	U.sig2_tot_CI95(y(3), :, :) = quantile(s2_tots, [0.025 0.975])';
      end

      if size(s2_tots, 2) == 1,
	s2s(:, j, :) = repmat(s2_tots, 1, 16);
      else
	s2s(:, j, :) = s2_tots;
      end
    end

    U.coeffs_sig2_tri(y(3), :) = glmfit(categs, s2s(:), 'normal');
  end

  V{x(3)} = U;
end

V = concat_structs(V);

%}}}

%
% 2. draw figure {{{
[~, sui] = unique(V.sig);

var_unk = cell(length(sui), 1);
var_exp = cell(length(sui), 1);

labels = cell(length(sui), 1);
sigs = cell(length(sui), 1);

for x = [sui [sui(2:end) - 1; slength(V)] (1:length(sui))']',
  i = x(1); j = x(2); k = x(3);

  Q = cumsum(V.coeffs_sig2_tri(i:j, :), 2);

  var_unk{k} = [sqrt(Q(:, 4)); 0];
  var_exp{k} = [max(0, diff(fliplr(sqrt(Q)), [], 2)); 0 0 0];

  labels{k} = mapacross(V.ch96(i:j), F.ch96, F.context);
  sigs{k} = V.sig{i};
end

%rename sigs to be human-readable
sigs = mapacross(sigs, {'APOBEC' 'ESO' 'LUNG' 'MSI' 'POLE' 'POLE_MSI' 'UV' 'VANILLA'}, ...
{'APOBEC' 'Eso.' 'Smoking' 'MSI' 'POLE' 'POLE+MSI' 'UV' 'CpG (Aging)'});
[sigs, si] = sort(sigs);
var_unk = var_unk(si);
var_exp = var_exp(si);
labels = labels(si);

figure(2); clf

ax_exp = gobjects(length(var_unk), 1);
ax_unk = gobjects(length(var_unk), 1);

cm = lines;

for i = 1:length(var_unk),
  %explained variance
  ax_exp(i) = axes;

  b = bar(var_exp{i}, 'stacked', 'BarWidth', 0.8);
  b(1).FaceColor = cm(1, :);
  b(2).FaceColor = cm(2, :);
  b(3).FaceColor = cm(3, :);

  ax_exp(i).XLim = [0.35 size(var_exp{i}, 1) - 0.35];
  ax_exp(i).YLim = [0 log(2)];

  ax_exp(i).YDir = 'reverse';

  ax_exp(i).XTickLabel = [];

  %unexplained variance
  ax_unk(i) = axes;

  b = bar(var_unk{i}, 'BarWidth', 0.8);
  b.FaceColor = 0.8*[1 1 1];

  ax_unk(i).XLim = [0.35 size(var_unk{i}, 1) - 0.35];
  ax_unk(i).YLim = [0 log(6)];

  ax_unk(i).XTickLabel = labels{i};

  ax_unk(i).XTickLabelRotation = 90;
  ax_unk(i).XAxis.FontSize = 10;
end

%get total width of figure
w = arrayfun(@(x) diff(x.XLim), ax_exp);
w_cs = cumsum([0; w(1:(end - 1))]);

%ticks
yt_exp = log(1:0.1:2);
yt_unk = log(1:0.2:6);

%reposition axes, apply formatting, add signature labels
for i = 1:length(var_unk),
  ax_exp(i).Position = [interval_remap(w_cs(i)/sum(w), 0, 1, 0.05, 0.95) + 0.005*(i - 1) 0.07 0.9*w(i)/sum(w) 0.2];
  ax_unk(i).Position = [interval_remap(w_cs(i)/sum(w), 0, 1, 0.05, 0.95) + 0.005*(i - 1) 0.45 0.9*w(i)/sum(w) 0.2*log(6)/log(2)];

  if i == 1,
    ax_exp(i).YLabel.String = {'exp(\sigma)', 'Explained'};
    ax_exp(i).YTickLabel = exp(yt_exp([1 4 6 8 11]));

    ax_unk(i).YLabel.String = {'exp(\sigma)', 'Unexplained'};
    ax_unk(i).YTickLabel = exp(yt_unk(1:5:end));
  else,
    ax_exp(i).YTickLabel = [];
    ax_unk(i).YTickLabel = [];
  end

  ax_exp(i).YTick = yt_exp([1 4 6 8 11]);
  ax_exp(i).YAxis.MinorTickValues = yt_exp;
  ax_exp(i).YAxis.MinorTick = 'on';
  ax_exp(i).YAxis.FontSize = 10;
  ax_exp(i).YGrid = 'on';
  ax_exp(i).YMinorGrid = 'on';

  ax_unk(i).YTick = yt_unk(1:5:end);
  ax_unk(i).YAxis.MinorTickValues = yt_unk;
  ax_unk(i).YAxis.MinorTick = 'on';
  ax_unk(i).YAxis.FontSize = 10;
  ax_unk(i).YGrid = 'on';
  ax_unk(i).YMinorGrid = 'on';

  ax_exp(i).XLabel.String = sigs{i};
  ax_exp(i).XLabel.Interpreter = 'none';
end

print('figures/variance_explained_NMF.svg', '-dsvg')

%postprocess SVG so that Illustrator can properly read fonts
% ssed -Ri "s/(.*font-family:')SansSerif('.*)/\1Helvetica\2/g" figures/variance_explained_NMF.svg

% }}}

%}}}

%
% Figure S5B {{{
clear

F = [];
F.file = direc('LNP_posteriors/signature_subcohorts/*/output/*.mat');
F = parsein(F, 'file', '.*/(.*?)/output/(\d+)-([ACGT]+)_(\d).*$', {'sig' 'c512' 'context' 'covars'});
F = makeapn(F);
F = reorder_struct(F, strcmp(F.sig, 'UV'));

F = sort_struct(F, {'c512' 'covars'});

[~, ui] = unique(F.c512);

var_exp_pct = NaN(length(ui), 3);
for x = [ui [ui(2:end) - 1; slength(F)] (1:length(ui))']',
  i = x(1); j = x(2); k = x(3);

  X1 = load(F.file{i});
  X2 = load(F.file{j});

  vep_raw = (1./X2.tau(1, 500:end) - 1./X1.tau(1, 500:end))./(1./X2.tau(1, 500:end));
  var_exp_pct(k, :) = [quantile(vep_raw, [0.1 0.9]) mean(vep_raw)];
end

[var_exp_pct, si] = sortrows(var_exp_pct, -3);

figure(2); clf
set(2, 'Position', [613 549 1310 420])
errorbar(1:size(var_exp_pct, 1), var_exp_pct(:, 3), diff(var_exp_pct(:, [1 3]), [], 2), diff(var_exp_pct(:, [3 2]), [], 2), 'linestyle', 'none', 'marker', '+')

ax = gca;
ax.YTick = -0.5:0.1:0.5;
ax.XTick = 1:size(var_exp_pct, 1);
ax.XLim = [0 size(var_exp_pct, 1) + 1];
ax.XTickLabel = F.context(ui(si));
ax.XTickLabelRotation = 45;
ax.XRuler.TickLength = [0.005 0.005];

line(xlim, [0 0], 'LineStyle', ':', 'Color', 'k', 'LineWidth', 2)

xlabel('Sequence context')
ylabel({'Fraction of additional variance' 'explained by XR-Seq covariate'})

print('figures/XR-seq.eps', '-depsc')

% }}}
