function [ax axT axO] = gene_alignment_plot(F, M, Mup, fnum, P)

demand_fields(F, {'gene' 'seq'})

%
%process params
if ~exist('P', 'var'), P = []; end
P = impose_default_value(P, 'text_highlight_index', 'cmat(:) >= 1');
P = impose_default_value(P, 'overlay_index', 'qmat < 0.5');
P = impose_default_value(P, 'overlay_color', 'ones(size(qmat))');
P = impose_default_value(P, 'x_axis_index', 1);
P = impose_default_value(P, 'xlim', []);

seqs = cat(1, F.seq{:});

%
%index coordinates WRT each gene
ngidx = seqs ~= '-';

coords = NaN(size(seqs));

for i = 1:slength(F),
  coords(i, ngidx(i, :)) = 1:nnz(ngidx(i, :));
end

%
%map mutations' coordinates
[~, pui, puj] = unique([M.mut.chr M.mut.pos M.mut.ch1536], 'rows');
ct = accumarray(puj, 1);

Mu = reorder_struct(M.mut, pui);
Mu.count = ct;
Mu.gene = M.gene.name(Mu.gene_idx);

Mu = multimapinto(Mu, Mup.mut, {'chr' 'pos' 'ch1536'}, {'chr' 'pos' 'ch1536'}, {'protein_change'});
Mu.protein_change(strcmp(Mu.protein_change, '')) = {'NA'};

Mu = reorder_struct(Mu, ismember(Mu.gene, F.gene));

Mu = parsein(Mu, 'protein_change', '([A-Z*])(\d+)([A-Z*])', {'aa_ref' 'aa_coord' 'aa_mut'});
Mu.aa_coord = str2double(Mu.aa_coord);

%
%matrix of q-values for each position
qmat = NaN(size(seqs));
for i = 1:slength(F),
  tmp = reorder_struct(Mu, strcmp(Mu.gene, F.gene{i}) & ~isnan(Mu.aa_coord));
  qv = accumarray(tmp.aa_coord, tmp.q(:, 3), [], @min, 1);
  qmat(i, :) = nansub(qv, coords(i, :), 1);
end

%
%matrix of counts for each position
cmat = NaN(size(seqs));
for i = 1:slength(F),
  tmp = reorder_struct(Mu, strcmp(Mu.gene, F.gene{i}) & ~isnan(Mu.aa_coord) & Mu.effect_idx ~= 5);
  ct = accumarray(tmp.aa_coord, tmp.count, [], @sum, 0);
  cmat(i, :) = nansub(ct, coords(i, :), 0);
end

%
%number of matches
eqidx = NaN([size(seqs) slength(F)]);
for i = 1:slength(F),
  eqidx(:, :, i) = bsxfun(@eq, seqs, seqs(i, :));
end
eqidx = sum(eqidx, 3);

%
%unique matches
uidx = NaN(size(seqs));
for i = 1:size(seqs, 2),
  [~, si] = sort(eqidx(:, i), 'descend');
  [~, ~, uidx(si, i)] = unique(seqs(si, i), 'stable');
end
uidx(eqidx == 1) = 0;

%
%draw plot {{{
[x, y] = meshgrid(1:size(seqs, 2), 1:size(seqs, 1));

figure(fnum); clf

%global xlim
if isempty(P.xlim), xl = [0 size(seqs, 2) + 1]; else, xl = P.xlim; end

%
%draw main plot {{{
mat = 4*(eqidx - 1) + uidx + 1;
image(mat)

colmap = [0.95 0.95 0.95; distinguishable_colors(length(unique(mat)) - 1, [1 1 1; 0 0 0; 1 0 1; 1 0 0; 1 0.67 0])];
colmap(end, :) = [1 0 1];
cmap = NaN(max(mat(:)), 3);
cmap(unique(mat(:)), :) = colmap;
colormap(cmap)

ax = gca;
ax.DataAspectRatio = [1 1 1];
ax.Position = [0.1 0 0.88 1];
ax.XTick = 1:size(eqidx, 2);
ax.XLim = [xl(1) + 0.5 xl(2) - 0.5];

%label x-axis WRT coordinates of specified alignment
n = P.x_axis_index;
ax.XTickLabel(ngidx(n, :)) = num2cellstr(coords(n, ngidx(n, :)));
ax.XTickLabel(~ngidx(n, :)) = {''};
ax.XTickLabel(ngidx(n, :) & mod(coords(n, :), 5) ~= 0) = {''};
ax.TickLength = [0 0];
ax.YTick = 1:slength(F);
ax.YTickLabel = F.gene;
ax.YAxis.FontSize = 10;

% }}}

%
%text axes {{{
axT = axes;
dummy = imagesc(NaN(size(qmat)));
dummy.Visible = 'off';

%index for *t*ext highlighting
tidx = eval(P.text_highlight_index);

%set x limits for text
lidx = x(:) < xl(2) & x(:) > xl(1);

text(x(~tidx & lidx), y(~tidx & lidx), seqs(~tidx & lidx), 'HorizontalAlignment', 'center', 'FontName', 'FixedWidth', 'FontSize', 8)
text(x(tidx & lidx), y(tidx & lidx), seqs(tidx & lidx), 'HorizontalAlignment', 'center', 'FontName', 'FixedWidth', 'FontSize', 8, 'Color', [1 0.67 0], 'FontWeight', 'bold')

axT.DataAspectRatio = ax.DataAspectRatio;
axT.Position = ax.Position;
axT.XTick = [];
axT.YTick = [];
axT.XLim = ax.XLim;
axT.Color = 'none';

% }}}

%
%overlay axes {{{
axO = axes;

%index for black *o*verlay
oidx = eval(P.overlay_index);

imagesc(eval(P.overlay_color), 'AlphaData', oidx)

axO.DataAspectRatio = ax.DataAspectRatio;
axO.Position = ax.Position;
axO.XTick = [];
axO.YTick = [];
axO.XLim = ax.XLim;
axO.Color = 'none';

colormap(axO, [0 0 0])

% }}}

% }}}
