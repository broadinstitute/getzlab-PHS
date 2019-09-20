function watercolor_plot(fignum, realdir, simdir)

%index data
F = [];
F.file = direc(realdir);
F = parsein(F, 'file', '.*/(.*?)/output/(\d+)-([ACGT]\(.->.\)[ACGT])_(\d).*$', {'sig' 'ch96' 'context' 'covars'});
F = makeapn(F);
F = reorder_struct(F, ~strcmp(F.sig, 'PANCAN') & F.covars == 3);

S = [];
[S.sig, sui] = unique(F.sig);

%marginals for each signature
U = cell(length(sui), 1);

%KDS confidence polygons
V = cell(length(sui), 1);

%2 std eigenvectors
W = cell(length(sui), 1);

pp = parpool(16);
for x = [sui [sui(2:end) - 1; slength(F)] (1:length(sui))']', 
  i = x(1); j = x(2); k = x(3);

  U{k} = cell(j - i + 1, 2);
  V{k} = cell(j - i + 1, 16);
  W{k} = cell(j - i + 1, 16);

  for l = i:j,
    %process real data (95CI regions and 2std bars)
    X = load(F.file{l});

    sig_tots = sqrt(bsxfun(@plus, 1./X.tau(1:16, 500:end), 1./X.tau(end, 500:end)))';
    mus = X.mu(:, 500:end)'; %XXX: account for covariates here too? this is mean sans covars.

    MG = cell(16, 1);
    CL = cell(16, 1);
    EV = cell(16, 1);
    parfor m = 1:16,
      xx = mus(:, m); yy = sig_tots(:, m);

      %marginal histogram
      MG{m} = {histc(xx, linspace(-16, -4, 500)) histc(yy, linspace(0, 2.5, 500))};

      %KDS smoothed confidence polygon
      [px py] = meshgrid(linspace(min(xx), max(xx), 100), linspace(min(yy), max(yy), 100));
      [f, xi] = ksdensity([xx yy], [px(:) py(:)]);

      C = contourc(px(1, 1:100), py(1:100, 1), reshape(f, 100, []), 40);

      cidx = 1;
      clidx = 1;
      clines = NaN;
      while true,
	if cidx > size(C, 2), break; end

	cidx0 = C(2, cidx);

	clines = C(:, 1 + cidx:(cidx + cidx0));

	mp = mean(inpolygon(xx, yy, clines(1, :), clines(2, :)));
	if mp >= 0.5 && mp <= 0.8, break; end

	cidx = 1 + cidx + cidx0;
	clidx = clidx + 1;
      end
      CL{m} = clines;

      %axes/center of 2 std confidence ellipse
      [ev lam] = eig(cov([xx yy]));
      evl = 2*ev*sqrt(lam);
      muh = mean([xx yy]);

      EV{m} = [evl; muh];
    end

    MG = reshape(cell2mat(cat(1, MG{:})), 500, []);
    U{k}(l - i + 1, :) = {sum(MG(:, 1:16), 2) sum(MG(:, 17:32), 2)};

    V{k}(l - i + 1, :) = CL(:);
    W{k}(l - i + 1, :) = EV(:);
  end
end
pp.delete

S.hist = U;
S.CR95 = V;
S.CV95 = W;
S.colors = [1 0 0; 0.74 0.68 0.052; 0 1 1; 0 1 0; 0 0 1; 0 0.5 0.5; 1 0.65 0; 0 0 0];
S.colors(S.colors == 1) = 1 - 2/255;
S.colors(S.colors == 0) = 2/255;
S.gfx = gobjects(slength(S), 1);

figure(fignum); clf
set(fignum, 'Position', [322 121 1278 726])
ax_plot = axes('Position', [0.05 0.3 0.7750 0.65]);
hold on

cof = dec2base(0:15, 3) - 48;
cof(cof == 1) = 2/255;
cof(cof == 2) = -2/255;

%plot 95CI ellipses
for s = 1:slength(S)
  for t = 1:size(S.CR95{s}, 1),
    for p = 1:16,
      x = S.CR95{s}{t, p}(1, :); y = S.CR95{s}{t, p}(2, :);
      parea = polyarea(x, y);
      pobj = patch(x, y, 0, ...
	   'FaceColor', S.colors(s, :) + cof(p, :), 'FaceAlpha', max(0.1, min(0.001/parea, 0.95)), ...
	   'EdgeColor', 'none');
      if t == 1 && p == 1, S.gfx(s) = pobj; end
    end
  end
end

%plot 2std lines
for s = 1:slength(S)
  for t = 1:size(S.CV95{s}, 1),
    for p = 1:16,
      q = S.CV95{s}{t, p};

      l1 = [q(3, :) - q(1:2, 1)'/2; ...
           q(3, :) + q(1:2, 1)'/2];
      l2 = [q(3, :) - q(1:2, 2)'/2; ...
           q(3, :) + q(1:2, 2)'/2];
      line([l1(:, 1); NaN; l2(:, 1)], [l1(:, 2); NaN; l2(:, 2)], 'Color', [S.colors(s, :) 0.2])
    end
  end
end

yt = log(0:12);

ax_plot.XLim = [-7*log(10) -2*log(10)];
ax_plot.YLim = [0 yt(end)];
ax_plot.YTick = yt;
ax_plot.YTickLabel = exp(yt);
ax_plot.XTick = [-7:-2]*log(10);
ax_plot.XMinorTick = 'on';
ax_plot.XAxis.MinorTickValues = log(logticks(1, 1e5)/1e7);
ax_plot.XTickLabel = 10.^([-7:-2] + 6);
ax_plot.XLabel.String = 'Median mutations per million sites [exp(\mu)]';
ax_plot.YLabel.String = 'Geometric standard deviation [exp(\sigma)]';
ax_plot.Box = 'on';
ax_plot.XGrid = 'on';
ax_plot.XMinorGrid = 'on';
ax_plot.YGrid = 'on';

%generate human-friendly legened
S.title = mapacross(S.sig, ...
{'APOBEC' ...
'ESO' ...
'LUNG' ...
'MSI' ...
'POLE' ...
'POLE_MSI' ...
'UV' ...
'VANILLA'}, ...
{'APOBEC' ...
'Esophageal' ...
'Smoking' ...
'MSI' ...
'POLE' ...
'POLE+MSI' ...
'UV' ...
'CpG (Aging)'});
[~, si] = sort(S.title);

leg = legend(S.gfx(si), S.title{si});
leg.Interpreter = 'none';

xl = ax_plot.XLim;
yl = ax_plot.YLim;

ax_mu = axes('Position', [0.05 0.01 0.775 0.2]);
hold on
ax_sig = axes('Position', [0.85 0.3 0.13 0.65]);
hold on
ax_sig.View = [90 90];
ax_sig.XDir = 'reverse';
for s = 1:slength(S),
  h = squeeze(sum(reshape(cell2mat(S.hist{s}), 500, [], 2), 2));
  axes(ax_mu)
  plot(linspace(xl(1), xl(2), 500), h(:, 1)./sum(h(:, 1)), 'Color', S.colors(s, :))
  axes(ax_sig)
  plot(linspace(yl(1), yl(2), 500), h(:, 2)./sum(h(:, 2)), 'Color', S.colors(s, :))
end
ax_mu.XLim = xl;
ax_mu.XTick = ax_plot.XTick;
ax_mu.XAxis.MinorTickValues = ax_plot.XAxis.MinorTickValues;
ax_mu.XMinorTick = 'on';
ax_mu.XAxisLocation = 'top';
ax_mu.XGrid = 'on';
ax_mu.XMinorGrid = 'on';
ax_mu.XTickLabels = [];
ax_mu.YTick = [];
ax_mu.Box = 'on';

ax_sig.YTick = [];
ax_sig.XTick = ax_plot.YTick;
ax_sig.XTickLabel = [];
ax_sig.XGrid = 'on';
ax_sig.Box = 'on';
