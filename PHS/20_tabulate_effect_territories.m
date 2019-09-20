%
% EFFECT ANALYSIS
%

%
% 1. tabulate effect territories for each pentamer, exomewide {{{

%look at overall effect territory for a given pentamer/base change combo, e.g.

%TCGTA -> (A/C/T) x (mis/non/syn)
%each pentamer is thus associated with nine values

clear

%load in target list
T = load_struct('ref/target_list.align75_filtered.txt');
T = makeapn(T);
T.tier = get_tiers(T.gene, 'ref/gencode/gene_tiers.v10.3.mat');
T.tier(isnan(T.tier)) = 4;

T = sort_struct(T, {'chr' 'start' 'end'});
[~, cui] = unique(T.chr);

%open C+E FWB
fwb = org.broadinstitute.cga.tools.seq.FixedWidthBinary('ref/gencode/c65e29/all.fwb');

%for each tier, tabulate effect vs. pentamer (non strand collapsed)
X = zeros([1024 29 6]);
for x = [cui [cui(2:end) - 1; slength(T)] (1:24)']'
  i = x(1); j = x(2); k = x(3);

  load(['ref/pentamer/chr' num2str(k) '.mat'])

  is_coding = false([size(categ, 1) 6]); 
  for l = i:j,
    is_coding(T.start(l):T.end(l), T.tier(l)) = true;
  end

  for t = 1:6,
    if nnz(is_coding(:, t)) == 0, continue; end %no genes in this tier in this chromosome

    ce = mod(double(fwb.get(T.chr(i), find(is_coding(:, t)))), 29) + 1;
    X(:, :, t) = X(:, :, t) + full(sparse(categ(is_coding(:, t)), ce, 1, 1024, 29));
  end

  fprintf('%d ', k);
end

%make collapse map
c1024 = dec2base(0:1023, 4) - 48;
c1024 = [(1:1024)' sum(bsxfun(@times, 3 - c1024(:, [1 3 2 5 4]), [4^4 4^3 4^2 4 1]), 2) + 1];
c512map = sparse(c1024(1:512, 2), ones(512, 1), c1024(1:512, 1));
c512map(1:512) = 1:512;
c512map = full(c512map);

efx = dec2base(0:26, 3) - 47;
efxrc = efx(:, [3 2 1]);
efxrcidx = sum(bsxfun(@times, efxrc - 1, [9 3 1]), 2) + 2;

%tabulate effect territories for each tier
effect_terrs = NaN([512 3 3 6]);
for t = 1:6,
  Xrc = X(513:end, :, t);
  Xrc(c512map(513:end), :) = Xrc;

  %reorder columns
  Xrc = Xrc(:, [1; efxrcidx; 29]);

  %add back
  X512 = X(1:512, :, t) + Xrc;

  %tabulate territories.  ch x effect x base
  for i = 1:512,
    idx = find(X512(i, 2:28));

    efxs = efx(idx, :);

    for j = 1:3,
      effect_terrs(i, :, j, t) = accumarray(efxs(:, j), X512(i, idx + 1)', [3 1]);
    end
  end
end

save('ref/pentamer_coding_effect_territories_by_tier.align75_filtered.mat', 'effect_terrs')

%}}}

%
% 2. tabulate effect territories for each trimer, for each gene {{{
clear

%use the MutSig gene definitions instead of the target list
load('ref/gencode/coverage_models.v3.mat')

X = reshape(C.gene_effect_cov(:, 3, :), C.ng, 29, []);
X(:, :, 65) = [];

%
%reverse complement 64 -> 32 {{{

c64 = dec2base(0:63, 4) - 48;
c64 = [(1:64)' sum(bsxfun(@times, 3 - c64(:, [1 3 2]), [4^2 4 1]), 2) + 1];
c32map = sparse(c64(1:32, 2), ones(32, 1), c64(1:32, 1));
c32map(1:32) = 1:32;
c32map = full(c32map);

%reorder columns
efx = dec2base(0:26, 3) - 47;
efxrc = efx(:, [3 2 1]);
efxrcidx = sum(bsxfun(@times, efxrc - 1, [9 3 1]), 2) + 1;

%add back
Xrc = X(:, :, 33:end);
Xrc(:, :, c32map(33:end)) = Xrc;
Xrc = Xrc(:, [efxrcidx; 28; 29], :);

X32 = X(:, :, 1:32) + Xrc;

%}}}

%tabulate effect territories for each gene (ch x effect x base x gene)
                           %ch e b g
effect_terrs_by_gene = NaN([32 4 3 C.ng]);
for g = 1:C.ng,
  for i = 1:32,
    %coding effects
    idx = find(X32(g, 1:27, i));

    efxs = efx(idx, :);

    for j = 1:3,
      effect_terrs_by_gene(i, 1:3, j, g) = accumarray(efxs(:, j), X32(g, idx, i)', [3 1]);
    end

    %splicing
    effect_terrs_by_gene(i, 4, :, g) = X32(g, 28, i);
  end
end

save('ref/trimer_coding_effect_territories_by_gene.mat', 'effect_terrs_by_gene')

%}}}
