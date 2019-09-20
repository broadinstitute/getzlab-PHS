%
%parse in runs {{{
clear

F = [];
F.results = direc('LNP_posteriors/signature_subcohorts_mutrate_split/*/output/*.mat');
F = parsein(F, 'results', '.*/([A-Z]+)_(hi|lo)/output/(\d+)-(.*)\.mat$', {'sig' 'rate' 'ch96' 'context'});
F = makeapn(F);
F = sort_struct(F, {'sig' 'ch96' 'rate'});
[~, ui] = unique_combos(F.sig, F.ch96);

H = rmfields(reorder_struct(F, ui), {'results' 'rate'});
H.results_hi = F.results(ui);
H.results_lo = F.results(ui + 1);

%}}}

%
% Figure S5A {{{

cmap = distinguishable_colors(16);

S = [];
S.sig = H.sig;
S.ch96 = H.ch96;
S.CR95 = cell(slength(H), 1);

%which contexts to include?
%eyeballed from plots of posterior vs. prior (see below); really, we should do this properly
%and compute some sort of mutual information (KL div, difference of entropy, JS div?)
S.include = false(slength(S), 16);

I = [];
I.sig = {'APOBEC' 'POLE' 'POLE' 'UV' 'UV'}';
I.ch96 = {96 64 95 93:95 96}';
I.pent = {9 1:16 1:16 1:16 [14 16]}';

for i = 1:slength(I),
  for j = find(strcmp(S.sig, I.sig{i}) & ismember(S.ch96, I.ch96{i}))',
    S.include(j, I.pent{i}) = true;
  end
end

%loop over contexts to include
for i = 1:slength(H),
  if ~any(S.include(i, :)), continue; end

  X_hi = load(H.results_hi{i});
  X_lo = load(H.results_lo{i});

  sig_tots_hi = sqrt(bsxfun(@plus, 1./X_hi.tau(1:16, 500:end), 1./X_hi.tau(end, 500:end)))';
  sig_tots_lo = sqrt(bsxfun(@plus, 1./X_lo.tau(1:16, 500:end), 1./X_lo.tau(end, 500:end)))';

  CL = cell(16, 1);

  %contour code excised from 50_signature_analysis
  for j = find(S.include(i, :)),
    xx = sig_tots_lo(:, j); yy = sig_tots_hi(:, j);

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
    CL{j} = clines;
  end

  S.CR95{i} = CL;
end

cols = [];
cols.sig = {'APOBEC' 'ESO' 'LUNG' 'MSI' 'POLE' 'POLE_MSI' 'UV' 'VANILLA'}';
cols.rgb = [1 0 0; 0.74 0.68 0.052; 0 1 1; 0 1 0; 0 0 1; 0 0.5 0.5; 1 0.65 0; 0 0 0];

S.colors = mapacross(S.sig, cols.sig, cols.rgb);

S.gfx = gobjects(slength(S), 1);

%scatterplot of joint sigma_hi vs. sigma_lo
figure(1); clf
hold on
for i = 1:slength(S),
  if ~any(S.include(i, :)), continue; end

  c = 0;
  for j = find(S.include(i, :)),
    x = S.CR95{i}{j}(1, :); y = S.CR95{i}{j}(2, :);
    parea = polyarea(x, y);

    %0.0043 is the minimum area of any joint
    %alph = 0.5*0.0043/parea;
    alph = max(0.1, min(0.25*0.0043/parea, 0.95));

    pobj = patch(x, y, 0, ...
	 'FaceColor', S.colors(i, :), 'FaceAlpha', alph, ...
	 'EdgeColor', 'none');

    if c == 0, S.gfx(i) = pobj; end
    c = c + 1;
  end
end
xlim([log(2) log(4.5)])
ylim([log(2) log(4.5)])

lt = 1:0.5:4.5;

ax = gca;

ax.XTick = log(lt);
ax.YTick = log(lt);
ax.XTickLabel = lt;
ax.YTickLabel = lt;
ax.Box = 'on';
line(xlim, xlim, 'LineStyle', '--', 'Color', 0.5*[1 1 1])
axis square

grid on

xlabel('exp(\sigma) non-hypermutants')
ylabel('exp(\sigma) hypermutants')

legend(S.gfx([4 8 13]), 'APOBEC', 'POLE', 'UV', 'Location', 'NorthWest')

print('figures/hypermut_sigmas.png', '-dpng', '-r600')

% }}}
