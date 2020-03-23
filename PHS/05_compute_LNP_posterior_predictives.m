%
% compute posterior predictives {{{

%
% APOBEC full covars vs. full+hairpin (Figure S5C) {{{

clear

%
% compute posterior predictives {{{

F = [];
F.file = direc('LNP_posteriors/signature_subcohorts/APOBEC/output/*.mat');
F = parsein(F, 'file', '.*/(.*?)/output/(\d+)-([ACGT]\(.->.\)[ACGT])_(\d).*$', {'sig' 'ch96' 'context' 'covars'});
F = makeapn(F);
F = reorder_struct(F, ismember(F.covars, [3 9]));

F = sort_struct(F, {'ch96' 'covars'});

PP = [];
[PP.ch96, ui] = unique(F.ch96);
PP.post_pred_pv = cell(slength(PP), 2);
PP.gene = cell(slength(PP), 1);
PP.protein_change = cell(slength(PP), 1);

%get protein coding effects
Mu = load2('tables/loci.prot_eff_harmonized.mat');
Mu = reorder_struct(Mu, ismember(Mu.ch96, F.ch96));

% load full MAF
% this time, we are considering all patients' mutations, not just those assigned to the APOBEC
% signature subcohort

M = loadM('mutation_data/MC3.align75.ICE_PoN-uniqued.M');

LuT = load_struct('ref/1536_LuT_v2.txt');
LuT = makeapn(LuT);

Lu = unique([LuT.c32 LuT.nbidx LuT.ch96], 'rows');

% isolate to just the relevant ch96's

for x = [ui [ui(2:end) - 1; slength(F)] (1:length(ui))']',
  i = x(1); j = x(2); k = x(3);

  ch96 = F.ch96(i);

  X1 = load(F.file{i});
  X2 = load(F.file{j}); %with hairpin

  % load mutations in this channel from PanCan MAF
  % get context 32 and newbase index corresponding to this ch96
  c32nb = Lu(Lu(:, 3) == ch96, 1:2);

  Mu = reorder_struct(M.mut, M.mut.c32 == c32nb(1));
  Mu.count = Mu.count_nb(:, c32nb(2));

  [MM C] = covar_annot(Mu, covars);
  [~, ~, C(:, end)] = unique(C(:, end));

  nzidx = MM > 0;
  

  %isolate to only sites with mutations
  %index will be same for both X1 and X2
  idx = X1.M > 0;

  X1.P.thin = 10;
  X1.P.pp_buffer_size = 260;
  X2.P.thin = 10;
  X2.P.pp_buffer_size = 260;

  [~, PP.post_pred_pv{k, 1}] = regr_post_pred(X1.mu, X1.tau, X1.Beta, X1.M(idx), X1.P.log_exposure, X1.C(idx, :), X1.P);
  [~, PP.post_pred_pv{k, 2}] = regr_post_pred(X2.mu, X2.tau, X2.Beta, X2.M(idx), X2.P.log_exposure, X2.C(idx, :), X2.P);

  %look up protein coding effects for these mutations
  Mu_ch = reorder_struct(Mu, Mu.ch96 == F.ch96(i));
  mmidx = multimap(X1.Mu, Mu_ch, {'chr' 'pos'});

  PP.protein_change{k} = nansub(Mu_ch.protein_change_harm2, mmidx, '');

  %genes for these mutations
  PP.gene{k} = X1.Mu.gene;
end

% %retroactively fix TBC1D12, which was missing hairpin covariates
% load('ref/apobec/all_TpC_sites.v2.0.mat', 'S')
% 
% i = 9; j = 10; k = 5;
% 
% X1 = load(F.file{i});
% X2 = load(F.file{j}); %with hairpin
% 
% %find TBC1D12 mutation
% idx = X1.M == 54;
% 
% look(X1.Mu, strcmp(X1.Mu.gene, 'TBC1D12'))
% 
% %          chr: 10
% %          pos: 96162370
% %      context: 526
% %          c65: 33
% %          c32: 32
% %         c512: 505
% %    is_coding: 1
% %         tier: 4
% %     gene_idx: 15997
% %      ref_idx: 3
% %     count_nb: 0 0 54
% %        count: 54
% %         gene: TBC1D12
% %            p: NaN
% 
% %fix covariates for that locus
% X2.C(idx, 5) = log10(X.site.relrate_exp(X.site.chr == 10 & X.site.pos == 96162370));
% 
% %compute posterior predictive at this site
% X2.P.thin = 10;
% X2.P.pp_buffer_size = 260;
% 
% idx_pp = X2.Mu.count == 54;
% 
% [~, PP.post_pred_pv{k, 2}(idx_pp, :, :)] = regr_post_pred(X2.mu, X2.tau, X2.Beta, X2.M(idx), X2.P.log_exposure, X2.C(idx, :), X2.P);

% retroactively add in RXRA, since it's not a T4 gene

% load in PanCan MAF
M = loadM('mutation_data/MC3.align75.ICE_PoN-uniqued.M');
hsidx = ismember(M.mut.gene_idx, find(ismember(M.gene.name, {'RXRA'}))) & M.mut.count > 10;

Muu = reorder_struct(M.mut, hsidx); Muu.context = Muu.c1025;

% get posteriors for this channel (94)
i = 9; j = 10; k = 5;

X1 = load(F.file{i});
X2 = load(F.file{j});

% add covariates for this locus
[MM C] = covar_annot(Muu, X2.P.covar_dir);
nzidx = MM > 0;
MM = MM(nzidx); C = C(nzidx, :);

% convert to categorical covariate
C(end) = find(unique(X1.Mu.context) == C(end));

[~, PP.post_pred_pv{k, 1}(end + 1, :, :)] = regr_post_pred(X1.mu, X1.tau, X1.Beta, MM, X1.P.log_exposure, C([1:4 6]), X1.P);
[~, PP.post_pred_pv{k, 2}(end + 1, :, :)] = regr_post_pred(X2.mu, X2.tau, X2.Beta, MM, X2.P.log_exposure, C, X2.P);
PP.gene{k}(end + 1) = {'RXRA'};

save('LNP_posteriors/signature_subcohorts/APOBEC/post_pred_hairpins.mat', 'PP')

% }}}

%
% draw Figure S5C {{{

% ! XXX ! following code will crash tmux

figure(4); clf
hold on
ax = gca;

cols = lines;

for i = 1:6,
  scatter(-log10(PP.post_pred_pv{i, 1}(:, 3, 1)), -log10(PP.post_pred_pv{i, 2}(:, 3, 1)), 'Marker', 'o', 'MarkerFaceColor', cols(1, :), 'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', 'none')
end

%q-cutoff
q_cut = -log10(4.88e-7);
patch([q_cut*[1 1] ax.XLim(2)*[1 1]], [ax.YLim fliplr(ax.YLim)], [0.5 0.5 0.5], 'EdgeColor', 'none', 'FaceAlpha', 0.25);
patch([ax.XLim fliplr(ax.XLim)], [q_cut*[1 1] ax.YLim(2)*[1 1]], [0.5 0.5 0.5], 'EdgeColor', 'none', 'FaceAlpha', 0.25);

%add text labels
pv = -log10(cell2mat(cellfun(@(x) x(:, 3, 1), PP.post_pred_pv, 'unif', 0)));

sigidx = sum(pv > q_cut, 2) > 0;

genes = cat(1, PP.gene{:});
pchange = cat(1, PP.protein_change{:});

pchange(sigidx & strcmp(genes, 'TBC1D12')) = {'M1*'};

textfit(pv(sigidx, 1) + 0.25, pv(sigidx, 2), strcat(genes(sigidx), {' '}, pchange(sigidx)))

%1:1 line
line(xlim, xlim, 'LineStyle', ':', 'Color', 'k')

%labels
ax.XTick = 0:2:14;
ax.YTick = 0:2:14;
ax.XTickLabel = ['1'; strcat('10^{-', ax.XTickLabel(2:end), '}')];
ax.YTickLabel = ['1'; strcat('10^{-', ax.YTickLabel(2:end), '}')];

ax.XLabel.String = 'p-value sans hairpin covariate';
ax.YLabel.String = 'p-value with hairpin covariate';

ax.Box = 'on';

print('figures/apo_scatter.svg', '-dsvg', '-painters')

% clean up SVG
% ssed -Ri -e "s/(.*font-family:')SansSerif('.*)/\1Helvetica\2/g" figures/apo_scatter.svg

% what's wrong with TBC1D12?  I thought this was one of the poster children of hairpin activity!
% is this issue that it's noncoding, so it didn't get any covariate annotations?
% if this is the case, we can re-annotate it manually and re-compute the posterior predictive just 
% for it -- I doubt its addition will severely alter the slope vector.

% and, while we're at it, check AHR

% for x = [ui [ui(2:end) - 1; slength(F)] (1:length(ui))']',
%   i = x(1); j = x(2); k = x(3);
% 
%   X1 = load(F.file{i});
%   X2 = load(F.file{j}); %with hairpin
% 
% %   %how many sites are affected by missing hairpin covariates, and what are their counts?
% %   count(X2.M(X2.C(:, 5) == 0))
% 
%   gidx = strcmp(X1.Mu.gene, 'AHR');
%   if ~any(gidx), continue; end
% 
%   idx = find(X1.M > 0);
% 
%   pidx = idx(gidx);
% 
%   for pp = pidx',
%     fprintf('%d: %0.2f %0.2f %0.2f %0.2f %0.2f   %d\n', k, X2.C(pp, 1:5), X2.M(pp, :))
%   end
% end
% 
% C = load('ref/cod_APOBEC.lawrence.hairpin_covars.align75_filtered/v1/chr10.mat')
% 
% C.covu(C.covar_mat(96162370, :)) %-Inf -> not in covariate track

%
% where is C3orf70? it doesn't appear in the APOBEC subcohort {{{

% M = loadM('mutation_data/MC3.align75.ICE_PoN_sig_annots.M');
% 
% signames = {'CpG' 'APOBEC' 'BRCA' 'smoking' 'ERCC2' 'MSI' 'UV' '8' 'AID/POLeta' 'POLE' 'alkylating' 'liver' 'APOBEC' '14' 'MSI' 'liver' 'eso' 'NB' 'PA' 'MSI' 'MSI' 'AA' '23' 'aflatoxin' 'hodgkin' 'MSI' 'KIRC' 'eso' '29' '30'};
% 
% c3idx = M.mut.chr == 3 & M.mut.pos == 184870595 & M.mut.ch96 == 95;
% 
% figure(6); clf
% imagesc(M.mut.sigmat(c3idx, :));
% ax = gca;
% ax.XTick = 1:30;
% ax.XTickLabel = signames;
% ax.XTickLabelRotation = 90;
% ax.CLim = [0 1];
% 
% ax.Title.String = 'C3orf70 S6L signature assignments';
% 
% cb = colorbar;
% cb.YLabel.String = 'Signature assignment probability';

% it doesn't appear because it's not assigned with probability >75% to APOBEC!

% }}}

% }}}

% }}}

%
% 2. pan-cancer {{{
clear

%load in runs
F = [];
F.file = direc('LNP_posteriors/pancancer/output/*.mat');
F = parsein(F, 'file', '.*/output/(\d+)-(.*)_\d\.mat$', {'ch' 'context'});
F = makeapn(F);

%loop over and extract p-values (for all loci)
pp = parpool(16);
for i = 1:slength(F),
  load(F.file{i});

  ints = min(length(M), ceil(linspace(1, length(M), 100)));
  ints = [ints(1:(end - 1)); [ints(2:(end - 1)) - 1 length(M)]]';
  
  P.thin = 10;
  P.pp_buffer_size = 260;

  post_pred = cell(size(ints, 1), 1);
  pvalues = cell(size(ints, 1), 1);
  parfor j = 1:size(ints, 1),
    r = ints(j, 1):ints(j, 2);
    [post_pred{j} pvalues{j}] = regr_post_pred(mu, tau, Beta, M(r), P.log_exposure, C(r, :), P);
  end

  post_pred = cat(1, post_pred{:});
  pvalues = cat(1, pvalues{:});

  post_pred = post_pred;
  pvalues = pvalues;
  Mu.p = pvalues(M > 0, :, 1);

  %overwrite previous version with updated posterior predictive/pvalues
  save(F.file{i}, 'M', 'P', 'Mu', 'C', 'Beta', 'mu', 'tau', 'logmarg_lik', 'post_pred', ...
                  'pvalues', 'logmarg_lik_unif')

  fprintf('%d ', i)
end 

%}}}

%}}}

%
% extract and save p-values from pan-cancer run
F = [];
F.file = direc(['LNP_posteriors/pancancer/output/*.mat']);
F = parsein(F, 'file', '.*/output/(\d+)-([ACGT]\(.->.\)[ACGT]).*\.mat$', {'ch96' 'context'});
F.ch96 = str2double(F.ch96);

[L P] = process_run(F);
save(['LNP_posteriors/pancancer/loci_pvalues.mat'], 'L');
