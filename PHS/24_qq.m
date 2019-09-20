%
% draw Q-Q plot (Figure S2C)


Q = [];
Q.file = direc('figures/p_mid*.mat');
Q.shortname = regexprep(Q.file, '.*_([A-Z]+)\.mat$', '$1');
Q.legendname = mapacross(Q.shortname, {'LNP' 'NB' 'UP' 'UWG'}, {'Log-normal-Poisson' 'Gamma-Poisson' 'Uniform Poisson' 'Uniform-within-gene'});
Q.color = mm_colormap;
Q.legobj = gobjects(slength(Q), 1);

figure(1); clf

%
%draw QQ {{{
ax_QQ = axes;
hold on

for i = 1:slength(Q),
  load(Q.file{i}, 'P', 'Z')

  %get relevant fieldname
  fn = fieldnames(P);
  f_rel = grep('^p.*mid$', fn);

  ps = sort(-log10([P.(f_rel{1}); Z.(f_rel{1})]));
  ps(ps > 16) = 16;
  np = length(ps);

  idx = ps > 2;
  n_cut = nnz(idx);

  Q.legobj(i) = scatter(-log10((n_cut:-1:1)/(np + 1)), ps(idx), 'Marker', '.', 'MarkerEdgeColor', Q.color(i, :));
  scatter(-log10((np:-10000:n_cut)/(np + 1)), ps(1:10000:(np - n_cut)), 'Marker', '.', 'MarkerEdgeColor', Q.color(i, :))
end

ax_QQ.Color = 'none';

% }}}

%
%draw CI's {{{
ax_CI = axes;

pts = round([np:(-np/1000):101  101:-1:1]);
beta95ci = NaN(length(pts), 2);
for y = [pts; 1:length(pts)],
  i = y(1); j = y(2);
  beta95ci(j, :) = -log10(icdf('beta', [0.025 0.975], i, np - i + 1));
end

x = as_column(-log10((np:-1:1)/(np + 1)));

fill([x(np - pts + 1); flipud(x(np - pts + 1))], [beta95ci(:, 2); flipud(beta95ci(:, 1))], [0.5 0.5 0.5], 'EdgeColor', 'none', 'FaceAlpha', 0.25)

% }}}

ax_CI = axes;
linkaxes([ax_QQ ax_CI])

xlim([0 8])
ylim([0 16])

leg2 = gobjects(2, 1);

%q < 0.1 region
leg2(1) = area(xlim, xlim + 1, ax_CI.YLim(2), 'LineStyle', ':', 'EdgeColor', 'k', 'FaceColor', [0.5 0.5 0.5], 'FaceAlpha', 0.25);

%1:1 line
leg2(2) = line(xlim, xlim, 'Color', 'k', 'LineStyle', '--');

ax_QQ.XAxis.Visible = 'off';
ax_QQ.YAxis.Visible = 'off';

uistack(ax_CI, 'bottom')

legend([Q.legobj; leg2], Q.legendname{:}, 'q < 0.1 cutoff', '1:1 line', 'Location', 'NorthWest')

ax_CI.XAxis.TickLabelFormat = '10^{-%g}';
ax_CI.YAxis.TickLabelFormat = '10^{-%g}';
ax_CI.XTickLabel{1} = '1';
ax_CI.YTickLabel{1} = '1';

xlabel('Uniform theoretical p-value quantiles')
ylabel('Observed p-value quantiles')

print('figures/QQ_for_publication.png', '-dpng', '-r500', '-painters')
