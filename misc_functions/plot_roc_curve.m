function ax = plot_roc_curve(G, pos_idx, neg_idx, X, q_range, fnum)

figure(fnum); clf
hold on

fprintf('FPR_____________ CI95___________________ TPR___ CI95_______ FDR___________ Algo\n')

for m = 1:slength(X),
  %
  %extract relevant q-values
  q = G.(X.fields{m});
  [~, si] = sort(q);
  pos_sort = pos_idx(si);
  neg_sort = neg_idx(si);

  %
  %ROC figure
  cs_ns = cumsum(neg_sort); n_tot_neg = nnz(neg_sort);
  x = log10(cs_ns/n_tot_neg);

  cs_ps = cumsum(pos_sort); n_tot_pos = nnz(pos_sort);
  y = cs_ps/n_tot_pos;

  X.roc_legends(m) = stairs(x, y, 'Color', X.colors(m, :), 'LineWidth', 0.5 + 0.75*(4 - m));
  qidx = NaN(length(q_range), 1);
  for z = [10.^q_range; 1:length(q_range)]
    i = z(1); j = z(2);
    qidx(j) = find(q(si) >= i, 1);
  end
  scatter(x(qidx), y(qidx), 'Marker', 'x', 'MarkerEdgeColor', X.colors(m, :))
  tx = strsplit(sprintf('10^{%d} ', q_range), ' '); tx(end) = [];
  text(x(qidx) + 0.05, y(qidx), tx, 'Color', X.colors(m, :))

  %print FPR/TPR (plus confidence intervals) and FDR at q-threshold of 0.1
  x2 = cumsum(neg_sort);
  a_fp = x2(qidx(2)); b_fp = nnz(neg_sort) - x2(qidx(2)) + 1;

  y2 = cumsum(pos_sort);
  a_tp = y2(qidx(2)); b_tp = nnz(pos_sort) - y2(qidx(2)) + 1;

  n_tot = a_tp + a_fp;

  fprintf('%0.4f (%0.3d/%0.3d) [%0.4e %0.4e] %0.4f [%0.2f %0.2f] %0.2f (%0.3d/%0.3d) %s\n', ...
          10.^x(qidx(2)), cs_ns(qidx(2)), n_tot_neg, icdf('beta', [0.025 0.975], a_fp, b_fp), y(qidx(2)), icdf('beta', [0.025 0.975], a_tp, b_tp), a_fp/n_tot, a_fp, n_tot, X.fields{m})
end

ax = gca;
ax.XTick = -3:0;
ax.XAxis.MinorTickValues = log10(logticks(0, 1e3)/1e3);
ax.XGrid = 'on';
ax.YGrid = 'on';
ax.XMinorGrid = 'on';
ax.XTickLabel = strsplit(sprintf('10^{%d} ', ax.XTick), ' ');
ax.Box = 'on';

%legend(X.roc_legends([2 3 4 1]), X.roc_legend_labels{:}, 'Location', 'SouthEast')
legend(X.roc_legends, X.roc_legend_labels{:}, 'Location', 'SouthEast')

xlabel('False Positive Rate')
ylabel('True Positive Rate')
