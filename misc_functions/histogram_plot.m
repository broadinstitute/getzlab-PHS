function [sp1 sp2] = histogram_plot(Xh, obs_counts, tot_sites, cmap, fnum)

oc = log10(obs_counts/sum(obs_counts));

%compute residuals (absolute error -- if we do multiple samples, this could also be RMS)
rs = cell(4, 1);
for m = 1:4,
  hs = sum(Xh{m});
  rs{m} = oc(1:length(hs)) - log10(hs./sum(hs));
end

len = min(cellfun(@length, rs));
rs_m = NaN(len, 4);

for m = 1:4,
  rs_m(:, m) = rs{m}(1:len);
end

figure(fnum); clf

%main plot
sp1 = subplot('Position', [0.1300 0.1100 0.7750 0.6]);
hold on

%area indicating fraction of <= 1 site
area([-0.5 999], -log10(tot_sites)*[1 1], 'BaseValue', -12, 'EdgeColor', 'none', 'FaceColor', 0.5*[1 1 1], 'FaceAlpha', 0.25)

%expected
for m = 1:4,
  hs = sum(Xh{m}); hs = hs./sum(hs);
  plot(0:(length(hs) - 1), log10(hs), 'Color', cmap(m, :))
  scatter(0:(length(hs) - 1), log10(hs), 'MarkerEdgeColor', cmap(m, :), 'Marker', 'o', 'MarkerFaceColor', cmap(m, :))
end

%observed
line([(0:999) - 0.4; (0:999) + 0.4; NaN(1, 1000)], ...
     [oc; oc; NaN(1, 1000)], ...
     'Color', 'k', 'LineWidth', 1.5, 'LineStyle', ':')


xlim([-0.5 25.5])

sp1.XTick = 0:25;

sp1.Box = 'on';
sp1.XLabel.String = 'Number of mutated patients';
sp1.YTickLabel = strsplit(sprintf('10^{%d} ', sp1.YTick), ' ');
sp1.YTickLabel{end - 1} = sprintf('%0.2f', 10^max(oc));
sp1.XGrid = 'on';
sp1.YGrid = 'on';

%residuals plot
sp2 = subplot('Position', [0.13 0.73 0.775 0.25]);
colormap(cmap)
b = bar(0:(len - 1), rs_m, 'EdgeColor', 'none', 'BarWidth', 1);
bl = b.BaseLine;
bl.LineStyle = ':';
xlim([-0.5 25.5])
sp2.XTick = sp1.XTick;
sp2.XTickLabels = [];
sp2.YTick = 0:3;
sp2.YTickLabel = strsplit(sprintf('10^{%d} ', sp2.YTick), ' ');
sp2.YAxis.MinorTickValues = log10(logticks(0, 1e3));
sp2.YMinorTick = 'on';
sp2.XGrid = 'on';
sp2.YGrid = 'on';
sp2.YMinorGrid = 'on';
sp2.YLabel.String = 'Residuals';

%link the x-axes of the two plots
linkaxes([sp1 sp2], 'x')
