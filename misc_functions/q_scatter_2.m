function ax = q_scatter_2(G, fieldx, fieldy)

demand_fields(G, {'bad_dNdS' 'true_pos'})

if ~exist('fieldx', 'var'), fieldx = 'unif'; end
if ~exist('fieldy', 'var'), fieldy = 'lnp'; end
min_q_field_x = ['min_q_' fieldx];
min_q_field_y = ['min_q_' fieldy];
max_sig_sil_field = ['max_sig_sil_' fieldx];

sig_idx = G.(min_q_field_x) <= 0.1 | G.(min_q_field_y) <= 0.1;

ax = axes;
hold on

%compute jitter amounts
xx = min(-log10(G.(min_q_field_x)), 8);
jidx = xx == 0;
xx(jidx) = xx(jidx) - 0.5*rand(nnz(jidx), 1);
jidx = xx == 8;
xx(jidx) = xx(jidx) + 0.5*rand(nnz(jidx), 1);

yy = min(-log10(G.(min_q_field_y)), 8);
jidx = yy == 0;
yy(jidx) = yy(jidx) - 0.5*rand(nnz(jidx), 1);
jidx = yy == 8;
yy(jidx) = yy(jidx) + 0.5*rand(nnz(jidx), 1);

%shadow genes that contain nothing significant
x = xx(~sig_idx);
y = yy(~sig_idx);

%downsample points for q = 1 
oidx = G.(min_q_field_x)(~sig_idx) == 1 & G.(min_q_field_y)(~sig_idx) == 1;
oidx_f = find(oidx);
oidx_ds = false(size(oidx)); oidx_ds(oidx_f(randsample(nnz(oidx), nnz(~oidx)))) = 1;

scatter(x(oidx_ds | ~oidx), y(oidx_ds | ~oidx), 'Marker', '.', 'MarkerEdgeColor', 0.5*[1 1 1])

%events that otherwise look OK
idx = sig_idx;
x = xx(idx);
y = yy(idx);
p_unk = scatter(x, y, 'Marker', '.', 'MarkerEdgeColor', 'k');

%genes with neutral dN/dS
idx = (grepm('^OR\d', G.gene) | G.bad_dNdS) & sig_idx;
x = xx(idx);
y = yy(idx);
p_dNdS = scatter(x, y, 'Marker', 's', 'MarkerEdgeColor', 'r');

%true positive genes
idx = G.true_pos;
x = xx(idx);
y = yy(idx);
p_truepos = scatter(x, y, 'Marker', 'o', 'MarkerEdgeColor', 'g');

pause(5)

ax.YRuler.Axle.LineStyle = 'none';
ax.XRuler.Axle.LineStyle = 'none';

xlim([-0.5 8.5])
ylim([-0.5 8.5])

line(xlim, [1 1], 'LineStyle', '--', 'Color', [0.6 0.6 0.6])
line([1 1], ylim, 'LineStyle', '--', 'Color', [0.6 0.6 0.6])

%
%break axes {{{

%x-breaks
line([-0.5 -0.1], -0.5*[1 1], 'Color', 'k', 'Clipping', 'off')
line([0.1 7.9], -0.5*[1 1], 'Color', 'k', 'Clipping', 'off')
line([8.1 8.5], -0.5*[1 1], 'Color', 'k', 'Clipping', 'off')

line([-0.1 -0.05 0.05 0.1], [-0.5 -0.3 -0.7 -0.5], 'Clipping', 'off', 'Color', 'k')
line([7.9 7.95 8.05 8.1], [-0.5 -0.3 -0.7 -0.5], 'Clipping', 'off', 'Color', 'k')

%y-breaks
line(-0.5*[1 1], [-0.5 -0.1], 'Color', 'k', 'Clipping', 'off')
line(-0.5*[1 1], [0.1 7.9], 'Color', 'k', 'Clipping', 'off')
line(-0.5*[1 1], [8.1 8.5], 'Color', 'k', 'Clipping', 'off')

line([-0.5 -0.3 -0.7 -0.5], [-0.1 -0.05 0.05 0.1], 'Clipping', 'off', 'Color', 'k')
line([-0.5 -0.3 -0.7 -0.5], [7.9 7.95 8.05 8.1], 'Clipping', 'off', 'Color', 'k')

% }}}

ax.XTick = [-0.25 1:7 8.25];
ax.YTick = [-0.25 1:7 8.25];
ax.XTickLabel = strsplit(sprintf('10^{-%d} ', 0:8), ' ');
ax.YTickLabel = strsplit(sprintf('10^{-%d} ', 0:8), ' ');
ax.XTickLabel{1} = 1;
ax.YTickLabel{1} = 1;
endlab = ['\leq' ax.XTickLabel{end - 1}];
ax.XTickLabel{end - 1} = endlab;
ax.YTickLabel{end - 1} = endlab;
ax.TickDir = 'out';

ax.PlotBoxAspectRatio = [1 1 1];
