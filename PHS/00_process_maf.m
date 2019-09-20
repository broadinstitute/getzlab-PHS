%
% add annotations necessary for downstream steps to MC3 calls
%

%
% 1. process raw MAF {{{

%
% NOTE: this requires the *controlled* MC3 MAF, which requires appropriate dbGaP authorization
% to access. Substituting "mc3.v0.2.8.CONTROLLED.maf" with a path to the non-controlled
% MAF should yield essentially equivalent results in downstream analyses.
%

% 1.1. take sSNVs called by MuTect
%      filter on all criteria except those that pertain to coding intervals

% awk -F'\t' 'BEGIN { OFS = FS } $111 ~ /MUTECT(\||$)/ && \
% ($109 == "NonExonic" || $109 == "NonExonic,bitgt" || $109 == "bitgt" || $109 == "PASS") \
% { print $5, $6, $10, $11, $13, $16, $41, $42 }' \
% mutation_data/mc3.v0.2.8.CONTROLLED.maf > \
% mutation_data/MC3_noheader.txt

% 1.2. take indels called by 2 or more callers
%      same filtering criteria as above.
%      note that we're only using these calls for dN/dS purposes (to infer TS genes)

% awk -F'\t' 'BEGIN { OFS = FS } $114 >= 2 && ($10 == "INS" || $10 == "DEL") && \
% ($109 == "NonExonic" || $109 == "NonExonic,bitgt" || $109 == "bitgt" || $109 == "PASS") \
% { print $5, $6, $9, $10, $11, $13, $16, $41, $42 }' \
% mutation_data/mc3.v0.2.8.CONTROLLED.maf > \
% mutation_data/MC3_indels_noheader.txt

%}}}

%
% 2. load and annotate MAF according to our transcript definitions {{{

%
% 2.1. sSNVs {{{

M = load_struct_noheader('mutation_data/MC3_noheader.txt', '%s%f%s%s%s%s%f%f', {'chr' 'pos' 'classification' 'ref' 'newbase' 'pat' 'n_ref' 'n_alt'});
M.chr = convert_chr(M.chr);

%
%get context/channels {{{
[~, ~, M.ch96, M.ch1536, M.c65, M.c1025] = get_channel_rates(M);

%RC c1025
c1024 = dec2base(0:1023, 4) - 48;
c1024 = [(1:1024)' sum(bsxfun(@times, 3 - c1024(:, [1 3 2 5 4]), [4^4 4^3 4^2 4 1]), 2) + 1];
c512map = sparse(c1024(1:512, 2), ones(512, 1), c1024(1:512, 1));

M.c512 = M.c1025;
idx = M.c512 > 512;
M.c512(idx) = c512map(M.c512(idx));

%RC c65
c64 = dec2base(0:63, 4) - 48;
c64 = [(1:64)' sum(bsxfun(@times, 3 - c64(:, [1 3 2]), [4^2 4 1]), 2) + 1]; 
c32map = sparse(c64(1:32, 2), ones(32, 1), c64(1:32, 1));

M.c32 = M.c65;
idx = M.c32 > 32;
M.c32(idx) = c32map(M.c32(idx));

%}}}

%
%annotate with genes/effects {{{
M = sort_struct(M, {'chr' 'pos'});

%map to targets
T = load_struct('ref/gencode/target_list.v1.txt');
T = makeapn(T);
T = sort_struct(T, {'chr' 'start' 'end'});

targ = map_mutations_to_targets_fast([M.chr M.pos], [T.chr T.start T.end]);
M.is_coding = targ > 0;
M.gene = repmat({''}, slength(M), 1);
M.gene(M.is_coding) = T.gene(targ(M.is_coding));

%add tiers
M.tier = get_tiers(M.gene, 'ref/gencode/gene_tiers.v10.3.mat');
M.tier(isnan(M.tier) & M.is_coding) = 4; 

%enforce effect 
M.effect = enforce_effect(M);

%}}}

%
%annotate with tumor types {{{
TS = load_struct('ref/TCGA/tissueSourceSite.tsv');
TS.StudyName = upper(TS.StudyName);
DS = load_struct('ref/TCGA/diseaseStudy.tsv');
DS.StudyName = upper(DS.StudyName);
TS = mapinto(TS, DS, 'StudyName', {'StudyAbbreviation'});
tcode = regexprep(M.pat, 'TCGA-(..)-.*', '$1');
M.ttype = mapacross(tcode, TS.TSSCode, TS.StudyAbbreviation);

%}}}

%
%convert to M struct and save {{{

%first, save results as full struct
save2(M, 'mutation_data/MC3.mafdir');

MM = [];
MM.mut = M;

%%pat
MM.pat = [];
[MM.pat.name, pui, MM.mut.pat_idx] = unique(MM.mut.pat);
MM.pat.nmut = accumarray(MM.mut.pat_idx, 1);
MM.pat.ttype = MM.mut.ttype(pui);

%%ttype
MM.ttype = [];
[MM.ttype.name, ~, MM.mut.ttype_idx] = unique(MM.mut.ttype);
MM.ttype.npat = full(sum(sparse(MM.mut.ttype_idx, MM.mut.pat_idx, 1) > 0, 2));
MM.ttype.nmut = accumarray(MM.mut.ttype_idx, 1);

%%gene
MM.gene = [];
[MM.gene.name, gui, MM.mut.gene_idx] = unique(MM.mut.gene);
MM.gene.npat = full(sum(sparse(MM.mut.gene_idx, MM.mut.pat_idx, 1) > 0, 2));
MM.gene.nmut = accumarray(MM.mut.gene_idx, 1);
MM.gene.tier = MM.mut.tier(gui);

%%allele
MM.allele = [];
[MM.allele.base, ~, tmp] = unique([MM.mut.ref; MM.mut.newbase]);
MM.mut.ref_idx = tmp(1:slength(MM.mut));
MM.mut.newbase_idx = tmp((slength(MM.mut) + 1):end);

%%effect
MM.effect = [];
[MM.effect.name, ~, MM.mut.effect_idx] = unique(MM.mut.effect);

%save
MM.mut = keep_fields(MM.mut, {'chr' 'pos' 'n_ref' 'n_alt' 'ch96' 'ch1536' 'c65' 'c32' 'c1025' 'c512' 'is_coding' 'tier' 'pat_idx' 'ttype_idx' 'gene_idx' 'ref_idx' 'newbase_idx' 'effect_idx'});
saveM(MM, 'mutation_data/MC3.M');

%}}}

%}}}

%
% 2.2. INDELs {{{

M = load_struct_noheader('mutation_data/MC3_indels_noheader.txt', '%s%f%s%s%s%s%s%f%f', {'chr' 'pos' 'effect' 'classification' 'ref' 'newbase' 'pat' 'n_ref' 'n_alt'});
M.chr = convert_chr(M.chr);

%remove mutations with bad annotations
M = reorder_struct(M, ~grepm('^(Mis|Non)', M.effect));

%
%annotate with genes {{{

M = sort_struct(M, {'chr' 'pos'});

%map to targets
T = load_struct('ref/gencode/target_list.v1.txt');
T = makeapn(T);
T = sort_struct(T, {'chr' 'start' 'end'});

targ = map_mutations_to_targets_fast([M.chr M.pos], [T.chr T.start T.end]);
M.is_coding = targ > 0;
M.gene = repmat({''}, slength(M), 1);
M.gene(M.is_coding) = T.gene(targ(M.is_coding));

%add tiers
M.tier = get_tiers(M.gene, 'ref/gencode/gene_tiers.v10.3.mat');
M.tier(isnan(M.tier) & M.is_coding) = 4; 

%}}}

%
%annotate with tumor types {{{
TS = load_struct('ref/TCGA/tissueSourceSite.tsv');
TS.StudyName = upper(TS.StudyName);
DS = load_struct('ref/TCGA/diseaseStudy.tsv');
DS.StudyName = upper(DS.StudyName);
TS = mapinto(TS, DS, 'StudyName', {'StudyAbbreviation'});
tcode = regexprep(M.pat, 'TCGA-(..)-.*', '$1');
M.ttype = mapacross(tcode, TS.TSSCode, TS.StudyAbbreviation);

%}}}

%
%convert to M struct and save {{{

MM = [];
MM.mut = M;

%%pat
MM.pat = [];
[MM.pat.name, pui, MM.mut.pat_idx] = unique(MM.mut.pat);
MM.pat.nmut = accumarray(MM.mut.pat_idx, 1);
MM.pat.ttype = MM.mut.ttype(pui);

%%ttype
MM.ttype = [];
[MM.ttype.name, ~, MM.mut.ttype_idx] = unique(MM.mut.ttype);
MM.ttype.npat = full(sum(sparse(MM.mut.ttype_idx, MM.mut.pat_idx, 1) > 0, 2));
MM.ttype.nmut = accumarray(MM.mut.ttype_idx, 1);

%%gene
MM.gene = [];
[MM.gene.name, gui, MM.mut.gene_idx] = unique(MM.mut.gene);
MM.gene.npat = full(sum(sparse(MM.mut.gene_idx, MM.mut.pat_idx, 1) > 0, 2));
MM.gene.nmut = accumarray(MM.mut.gene_idx, 1);
MM.gene.tier = MM.mut.tier(gui);

%%allele
MM.allele = [];
[MM.allele.base, ~, tmp] = unique([MM.mut.ref; MM.mut.newbase]);
MM.mut.ref_idx = tmp(1:slength(MM.mut));
MM.mut.newbase_idx = tmp((slength(MM.mut) + 1):end);

%%effect
MM.effect = [];
[MM.effect.name, ~, MM.mut.effect_idx] = unique(MM.mut.effect);

%save
MM.mut = keep_fields(MM.mut, {'chr' 'pos' 'n_ref' 'n_alt' 'is_coding' 'tier' 'pat_idx' 'ttype_idx' 'gene_idx' 'ref_idx' 'newbase_idx' 'effect_idx'});
saveM(MM, 'mutation_data/MC3_indels.M');

%}}}

%
%concatenate with sSNVs, convert to M struct, and save {{{
MS = load2('mutation_data/MC3.mafdir');

MM = [];
MM.mut = concat_structs_keep_all_fields({M MS});

%%pat
MM.pat = [];
[MM.pat.name, pui, MM.mut.pat_idx] = unique(MM.mut.pat);
MM.pat.nmut = accumarray(MM.mut.pat_idx, 1);
MM.pat.ttype = MM.mut.ttype(pui);

%%ttype
MM.ttype = [];
[MM.ttype.name, ~, MM.mut.ttype_idx] = unique(MM.mut.ttype);
MM.ttype.npat = full(sum(sparse(MM.mut.ttype_idx, MM.mut.pat_idx, 1) > 0, 2));
MM.ttype.nmut = accumarray(MM.mut.ttype_idx, 1);

%%gene
MM.gene = [];
[MM.gene.name, gui, MM.mut.gene_idx] = unique(MM.mut.gene);
MM.gene.npat = full(sum(sparse(MM.mut.gene_idx, MM.mut.pat_idx, 1) > 0, 2));
MM.gene.nmut = accumarray(MM.mut.gene_idx, 1);
MM.gene.tier = MM.mut.tier(gui);

%%allele
MM.allele = [];
[MM.allele.base, ~, tmp] = unique([MM.mut.ref; MM.mut.newbase]);
MM.mut.ref_idx = tmp(1:slength(MM.mut));
MM.mut.newbase_idx = tmp((slength(MM.mut) + 1):end);

%%effect
MM.effect = [];
[MM.effect.name, ~, MM.mut.effect_idx] = unique(MM.mut.effect);

%save
MM.mut = keep_fields(MM.mut, {'chr' 'pos' 'n_ref' 'n_alt' 'ch96' 'ch1536' 'c65' 'c32' 'c1025' 'c512' 'is_coding' 'tier' 'pat_idx' 'ttype_idx' 'gene_idx' 'ref_idx' 'newbase_idx' 'effect_idx'});
saveM(MM, 'mutation_data/MC3_sSNVs-indels.M');

%}}}

%}}}

%}}}

%
% 3. restrict to unique regions (align75).  apply ICE PoN. {{{
clear

M = loadM('mutation_data/MC3.M');
M.mut = sort_struct(M.mut, {'chr' 'pos'});

T = load_struct('ref/target_list.align75_filtered.txt');
T = makeapn(T);
T = sort_struct(T, {'chr' 'start' 'end'});

targ = map_mutations_to_targets_fast([M.mut.chr M.mut.pos], [T.chr T.start T.end]);
M.mut = reorder_struct(M.mut, targ > 0);

%
%apply ICE PoN (need appropriate dbGaP access to download)
M.mut.pon = get_pon(M.mut.chr, M.mut.pos, 'ref/pon/ice_pon.bin') + get_pon(M.mut.chr, M.mut.pos, 'ref/pon/ice_pon2.bin'); 
M.mut.pon_ll = get_loglikelihood_from_pon_vars(pon, M.mut.n_alt, M.mut.n_ref);

idx = M.mut.pon_ll > -2.5;

%
%what sites get removed? {{{
Mu = reorder_struct(M.mut, idx);
[~, pui, puj] = unique([Mu.chr Mu.pos], 'rows');
Mu = reorder_struct(Mu, pui);
h = accumarray(puj, 1);
Mu.count = h;
Mu = sort_struct(Mu, 'count', -1);
Mu.gene = M.gene.name(Mu.gene_idx);
pr(Mu, {'chr' 'pos' 'n_ref' 'n_alt' 'effect_idx' 'count' 'gene' 'pon' 'pon_ll'}, 1:20)

% }}}

M.mut = reorder_struct(M.mut, ~idx);

%
%save
saveM(M, 'mutation_data/MC3.align75.ICE_PoN.M');

% }}}

%
% 4. make uniqued MAFS {{{

%unique on both positions and changes

%
% 3.1. PanCan {{{
clear

M = loadM('mutation_data/MC3.M');

M.mut.newbase_idx_rc = M.mut.newbase_idx;
idx = M.mut.c1025 > 512;
M.mut.newbase_idx_rc(idx) = 5 - M.mut.newbase_idx(idx); 
ch_base = dec2base(M.mut.c512 - 1, 4, 5) - 47; M.mut.ch_base = ch_base(:, 1);

nbidx = [0 1 2 3; ...
         1 0 2 3; ...
         1 2 0 3; ...
         1 2 3 0];

[~, pui, puj] = unique([M.mut.chr M.mut.pos], 'rows');
h = full(sparse(puj, nbidx(sub2ind(size(nbidx), M.mut.ch_base, M.mut.newbase_idx_rc)), 1));
M.mut = keep_fields(reorder_struct(M.mut, pui), {'chr' 'pos' 'c1025' 'c65' 'c32' 'c512' 'is_coding' 'tier' 'gene_idx' 'ref_idx'});
M.mut.count_nb = h;
M.mut.count = sum(M.mut.count_nb, 2);

%save
saveM(M, 'mutation_data/MC3-uniqued.M')

%}}}

%
% 3.11. PanCan, align75 + ICE PoN filtered {{{
clear

M = loadM('mutation_data/MC3.align75.ICE_PoN.M');

M.mut.newbase_idx_rc = M.mut.newbase_idx;
idx = M.mut.c1025 > 512;
M.mut.newbase_idx_rc(idx) = 5 - M.mut.newbase_idx(idx); 
ch_base = dec2base(M.mut.c512 - 1, 4, 5) - 47; M.mut.ch_base = ch_base(:, 1);

nbidx = [0 1 2 3; ...
         1 0 2 3; ...
         1 2 0 3; ...
         1 2 3 0];

[~, pui, puj] = unique([M.mut.chr M.mut.pos], 'rows');
h = full(sparse(puj, nbidx(sub2ind(size(nbidx), M.mut.ch_base, M.mut.newbase_idx_rc)), 1));
M.mut = keep_fields(reorder_struct(M.mut, pui), {'chr' 'pos' 'c1025' 'c65' 'c32' 'c512' 'is_coding' 'tier' 'gene_idx' 'ref_idx'});
M.mut.count_nb = h;
M.mut.count = sum(M.mut.count_nb, 2);

%save
saveM(M, 'mutation_data/MC3.align75.ICE_PoN-uniqued.M')

% }}}

%
% 3.31. Signature subcohorts {{{
clear

M = loadM('mutation_data/MC3.align75.ICE_PoN.M');

M.mut.newbase_idx_rc = M.mut.newbase_idx;
idx = M.mut.c1025 > 512;
M.mut.newbase_idx_rc(idx) = 5 - M.mut.newbase_idx(idx); 
ch_base = dec2base(M.mut.c512 - 1, 4, 5) - 47; M.mut.ch_base = ch_base(:, 1);

nbidx = [0 1 2 3; ...
         1 0 2 3; ...
         1 2 0 3; ...
         1 2 3 0];

%
%Nick's definition of POLE+MSI {{{

% from Haradhvala et al., "Distinct mutational signatures characterize concurrent loss of polymerase proofreading and mismatch repair", Nature Communications 2018, DOI https://doi.org/10.1038/s41467-018-04002-4

C = load_struct('ref/TCGA/UCEC.signature_based_classifications.txt');

pidx = strcmp(mapacross(M.pat.name, C.name, C.classification), 'POLE+MSI');
midx = ismember(M.mut.pat_idx, find(pidx));

Ms = M; 
Ms.mut = reorder_struct(Ms.mut, midx); 
Ms.npat = nnz(pidx);

Ms.pat = keep_fields(M.pat, {'name' 'ttype'});

%draw lego plot
Ms.mut.context1025 = Ms.mut.c1025;
Ms.mut.context65 = Ms.mut.c65;

%save non-uniqued MAF
saveM(Ms, 'mutation_data/MC3.align75.ICE_PoN.sig_cohorts/POLE_MSI.M')

figure(1); clf
set(1, 'Pos', [520 678 1120 420])
subplot('Position', [0.05 0.05 0.45 0.95])
lego5(Ms.mut)
zl = zlim;
text(1, 1, zl(2), 'POLE+MSI')
subplot('Position', [0.5 0.05 0.45 0.95])
lego(Ms.mut)

print('figures/signature_cohort_legos/POLE_MSI.png', '-dpng')

%unique MAF
[~, pui, puj] = unique([Ms.mut.chr Ms.mut.pos], 'rows');
h = full(sparse(puj, nbidx(sub2ind(size(nbidx), Ms.mut.ch_base, Ms.mut.newbase_idx_rc)), 1));
Ms.mut = keep_fields(reorder_struct(Ms.mut, pui), {'chr' 'pos' 'c1025' 'c65' 'c32' 'c512' 'is_coding' 'tier' 'gene_idx' 'ref_idx'});
Ms.mut.count_nb = h;
Ms.mut.count = sum(Ms.mut.count_nb, 2);

%save
saveM(Ms, 'mutation_data/MC3.align75.ICE_PoN-uniqued.sig_cohorts/POLE_MSI.M')

%}}}

%
% 3.31.2. NMF {{{

%
% 0. R code to import Bayesian NMF signature assignments from Jaegil Kim {{{

% setwd("mutation_data")
% load("mutation_data/maf.prob.pan_noDNP.SNV.projection_with_local.RData")
% 
% maf = as.data.frame(maf);
% 
% maflite = maf[, is.element(colnames(maf), c("Chromosome", "Start_Position", "patient_name", "context96.num")) | grepl("^Sig", colnames(maf))];
% 
% write.table(maflite, file="mutation_data/MC3_Jaegil_sig_annots.maf", row.names = FALSE, sep = '\t', quote = FALSE)

% }}}

%
% 1. load in MAF with signature assignments; map to our MAF {{{

Ms = load_struct('mutation_data/MC3_Jaegil_sig_annots.maf');
Ms = makeapn(Ms);
Ms.Chromosome = convert_chr(Ms.Chromosome);

Ms.pat_idx = listmap(Ms.patient_name, M.pat.name);

Ms.sigmat = struct2matrix(keep_fields(Ms, grep('^Signature\d+$', fieldnames(Ms))));

% map signature assignments into our callset

M.mut = multimapinto(M.mut, Ms, {'chr' 'pos' 'pat_idx' 'ch96'}, {'Chromosome' 'Start_Position' 'pat_idx' 'context96num'}, {'sigmat'});

M.mut.pat_ch_idx = (M.mut.pat_idx - 1)*96 + M.mut.ch96;
M.mut = sort_struct(M.mut, 'pat_ch_idx');

M.pat.sigmat = NaN(slength(M.pat), 96, 30);

[pu, pui] = unique(M.mut.pat_ch_idx);

for x = [pui [pui(2:end) - 1; slength(M.mut)] pu]',
  i = x(1); j = x(2); k = x(3);

  ch96 = mod(k - 1, 96) + 1;
  pat_idx = ceil(k/96);

  smat = unique(M.mut.sigmat(i:j, :), 'rows');
  smat(isnan(smat(:, 1)), :) = [];

  if size(smat, 1) == 0, continue; end

  M.pat.sigmat(pat_idx, ch96, :) = smat;
end

%save

saveM(M, 'mutation_data/MC3.align75.ICE_PoN_sig_annots.M');

% }}}

%
% 3. context/patient combos with high probabilistic assignments {{{

% in patient p, if relevant channel c is assigned to corresponding signature with prob >= 75%,
% then include all mutations of c in p

clear

M = loadM('mutation_data/MC3.align75.ICE_PoN_sig_annots.M');

X = [];
X.names = {'VANILLA' 'UV' 'LUNG' 'ESO' 'APOBEC' 'POLE' 'MSI'}';
X.summary = [0 2.^[0 1 3:6]]';

Y = generate_signature_subcohort_definitions('mutation_data/MC3.align75.ICE_PoN-uniqued.sig_cohorts/*.M');

%add back POLE-MSI manually
X.names{8} = 'POLE_MSI';
X.summary(8) = NaN;

X = mapinto(X, Y, 'names');

%add COSMIC annotations
X.cosmic = cell(slength(X), 1);
X.cosmic(1) = {1};
X.cosmic(2) = {7};
X.cosmic(3) = {4};
X.cosmic(4) = {17};
X.cosmic(5) = {[2 13]};
X.cosmic(6) = {10};
X.cosmic(7) = {6};
X.cosmic(8) = {14};

%some stuff for plotting
nbidx = [0 1 2 3; ...
         1 0 2 3; ...
         1 2 0 3; ...
         1 2 3 0];

%for plotting "jaEGO" plots
lego_colors = [1 1 0;0 0.7 0.7;1 0 0;0.1 0.8 0.1;0 0.2 0.8;0.5 0.3 0.7;];
lego_colors = lego_colors([5 4 6 2 3 1], :);

%for 3D LEGO plots
lego_colors_3d = get_LEGO_colors;

bc = {'A->C' 'A->G' 'A->T' 'C->A' 'C->G' 'C->T'};
B = {'A' 'C' 'G' 'T'};

tmp = dec2base((1:96) - 1, 4) - 48;
name96 = strcat(B(tmp(:, end - 1) + 1), '(', bc(sum(bsxfun(@times, tmp(:, 1:end - 2), [4 1]), 2) + 1), ')', B(tmp(:, end) + 1))';

%rename name96 to convert to Mike's LEGO plot definitions
name96_mike = cellfun(@(x) [x(3) ' in ' x(1) '_' x(end) ' ->' x(6)], name96, 'unif', 0);

%patients x channels counts
pat_ch96 = full(sparse(M.mut.pat_idx, M.mut.ch96, 1));

% Supp. Figure 4
for i = 1:slength(X),
  %
  %get relevant channels for this signature {{{
  ch96 = cell(length(X.output{i}), 1);

  for j = 1:length(X.output{i}),
    Y = [];
    Y.c512 = X.output{i}{j}(:, 1); 
    b = dec2base(Y.c512 - 1, 4, 5) - 48;

    Y.nb = X.output{i}{j}(:, 2);

    ch96{j} = unique(1 + [b(:, [2 3]) Y.nb - 1 b(:, 1)]*[1 4 16 48]');
  end
  ch96 = cat(1, ch96{:});

  % }}}

  %
  % index patients/channels strongly assigned to this signature {{{

  patmat = sum(M.pat.sigmat(:, ch96, X.cosmic{i}) > 0.75, 3);
  
  %to display all assigned channels (not just ones we've preselected):
  %patmat = sum(M.pat.sigmat(:, :, X.cosmic{i}) > 0.75, 3);

  [a b] = find(patmat);

  % }}}

  %
  %save MAFs {{{

  %note that in this case, we are not saving all mutations in relevant patients, only mutations
  %in relevant channels

  midx = ismember([M.mut.pat_idx M.mut.ch96], [a ch96(b)], 'rows');

  Ms = M;
  Ms.mut = reorder_struct(Ms.mut, midx);
  Ms.npat = sum(sparse(ch96(b), a, 1), 2);

  %unique
  [~, pui, puj] = unique([Ms.mut.chr Ms.mut.pos], 'rows');
  h = full(sparse(puj, nbidx(sub2ind(size(nbidx), Ms.mut.ch_base, Ms.mut.newbase_idx_rc)), 1));
  Ms.mut = keep_fields(reorder_struct(Ms.mut, pui), {'chr' 'pos' 'c1025' 'c65' 'c32' 'c512' 'is_coding' 'tier' 'gene_idx' 'ref_idx'});
  Ms.mut.count_nb = h;
  Ms.mut.count = sum(Ms.mut.count_nb, 2);

  %save
  saveM(Ms, ['mutation_data/MC3.align75.ICE_PoN-uniqued.sig_cohorts_NMF/' X.names{i} '.M'])

  % }}}

  %
  %save spectrum plots {{{

  %for supp. figure: LEGO/jaEGO plots of all the patients with at least one assignment; distinguish
  %fraction of counts assigned vs. total counts

  h = histc(M.mut.ch96(ismember(M.mut.pat_idx, a)), 1:96);
  h_a = accumarray(ch96(b), pat_ch96(sub2ind(size(pat_ch96), a, ch96(b))), [96 1]);

  %if we are displaying all assigned channels
  %h_a = accumarray(b, pat_ch96(sub2ind(size(pat_ch96), a, b)), [96 1]);

  ax = gobjects(6, 1);

  %jaEGO

  figure(i*100); clf
  set(gcf, 'Position', [439 859 1121 213])
  for x = [1:16:96; 16:16:96; 1:6],
    u = x(1); v = x(2); w = x(3);

    ax(w) = subplot('Position', [0.1 + 0.85*((w - 1)/6) 0.1 0.8/6 0.8]);
    hold on
    bar(h(u:v), 'FaceColor', lego_colors(w, :), 'FaceAlpha', 0.5, 'EdgeColor', 'none')
    bar(h_a(u:v), 'FaceColor', lego_colors(w, :))

    title(name96{u}(3:6))

    if w > 1, ax(w).YTickLabel = []; end
  end

  ax(1).YLim = [0 1.1*max(h)];
  ax(1).XLim = [0 17]; 
  linkaxes(ax)

  %LEGO

  %map to get ch96 in correct order for bar3
  h_mat = factor2lego(h', name96_mike);
  h_a_mat = factor2lego(h_a', name96_mike);

  %index zero assigned mutation counts
  [a, b] = find(h_a_mat == 0);
  h_a_zidx = full(sparse(a, b, 1));

  % Figure S4
  figure(i*1000); clf
  b_upper = bar3_with_colors(h_mat - h_a_mat, lego_colors_3d);
  hold on
  b_lower = bar3_with_colors(h_a_mat, lego_colors_3d);

  %increment z coordinates of upper bars; make upper bars transparent
  for r = 1:12,
    for c = (1:6:48) - 1,
      zidx = sub2ind([48 4], [1 1 2 2 3 3 4 4 5 5] + c, [2 3 1 4 1 4 2 3 2 3]);

      l_height = b_lower(r).ZData(c + 2, 2);
      u_height = b_upper(r).ZData(c + 2, 2);
      
      b_upper(r).ZData(zidx) = l_height;
      b_upper(r).ZData((c + 2):(c + 3), 2:3) = u_height + l_height;

      %remove lower bars with zero mutations
      if h_a_zidx(c/6 + 1, r),
	b_lower(r).ZData((c + 1):(c + 6), :) = NaN;
      end

      %remove occluded edges of upper bars to increase legibility
      bidx = sub2ind([48 4], [1 5 2] + c, [3 3 4]);
      b_upper(r).ZData(bidx) = NaN;
    end

    b_upper(r).FaceAlpha = 0.2;
    b_upper(r).EdgeAlpha = 0.2;
    b_upper(r).BackFaceLighting = 'unlit';
  end

  ax = gca;
  zmax = ceil(max(h)/1e4)*1e4; 
  ax.ZLim = [0 zmax];
  ax.XTick = [];
  ax.YTick = [];
  %ax.ZTick = 0:1e4:zmax;
  ax.ZAxis.Exponent = 4;
  ax.XAxis.Visible = 'off';
  ax.YAxis.Visible = 'off';

  ax.ZLabel.String = 'No. of mutations';
  ax.ZLabel.Rotation = 90;

  ax.Title.String = X.title{i};
  ax.Title.Interpreter = 'none';

  print(['figures/signature_cohort_legos_NMF/' X.names{i} '.svg'], '-dsvg')

  % }}}
end

%have to do some postprocessing to the SVGs so that Illustrator can read them

% find figures/signature_cohort_legos_NMF -name "*.svg" | while read -r i; do
% ssed -Ri -e "s/(.*font-family:')SansSerif('.*)/\1Helvetica\2/g" -e "s/(.*font-family:')mwb_cmsy10('.*)/\1Arial Unicode\2/g" -e "s/>#</>\&#x00D7;</g" $i
% done

% }}}

%}}}

%}}}

%
% 3.32. Signature subcohorts, split by mutation rate {{{

clear

M = loadM('mutation_data/MC3.align75.ICE_PoN_sig_annots.M');

%generate signature definitions
X = generate_signature_subcohort_definitions('mutation_data/MC3.align75.ICE_PoN-uniqued.sig_cohorts_NMF/*.M');

%add COSMIC annotations
X.cosmic = cell(slength(X), 1);
X.cosmic(8) = {1};
X.cosmic(7) = {7};
X.cosmic(3) = {4};
X.cosmic(2) = {17};
X.cosmic(1) = {[2 13]};
X.cosmic(5) = {10};
X.cosmic(4) = {6};
X.cosmic(6) = {14};

pat_ch96 = full(sparse(M.mut.pat_idx, M.mut.ch96, 1));

%
%load territories (for tabulating mutation rates)
load('ref/gene_list.align75_filtered.territories.mat', 'G')

%lookup tables for ch96 -> c32 x 3
tmp = dec2base(0:95, 4) - 48;
c96_32map = [tmp(:, 1:2)*[4 1]' > 2 tmp(:, 3:4)]*[16 4 1]' + 1;
%c96_bc = mod(ceil((1:96)/16) - 1, 3)' + 1;

%total number of trimers in exome
N32 = sum(G.terr32);
N96 = N32(c96_32map);

nbidx = [0 1 2 3; ...
         1 0 2 3; ...
         1 2 0 3; ...
         1 2 3 0];

%
% generate split MAFs for each signature {{{

for i = 1:slength(X),
  %
  %get relevant channels for this signature {{{
  ch96 = cell(length(X.output{i}), 1);

  for j = 1:length(X.output{i}),
    Y = [];
    Y.c512 = X.output{i}{j}(:, 1); 
    b = dec2base(Y.c512 - 1, 4, 5) - 48;

    Y.nb = X.output{i}{j}(:, 2);

    ch96{j} = unique(1 + [b(:, [2 3]) Y.nb - 1 b(:, 1)]*[1 4 16 48]');
  end
  ch96 = cat(1, ch96{:});

  % }}}

  %
  % index patients/channels strongly assigned to this signature

  patmat = sum(M.pat.sigmat(:, ch96, X.cosmic{i}) > 0.75, 3);

  [a b] = find(patmat);

  %
  % generate set of non-uniqued mutations

  midx = ismember([M.mut.pat_idx M.mut.ch96], [a ch96(b)], 'rows');

  Ms = M;
  Ms.mut = reorder_struct(Ms.mut, midx);

  %
  % split by patients' mutation rate {{{

  % which channels to include in rate calculation will depend on the patient
  pat_ch = accumarray(a, ch96(b), [9023 1], @(x) {x}); 

  up = unique(a);

  patrates = NaN(length(up), 1);
  patcounts = NaN(length(up), 1);

  for p = [up'; 1:length(up)],
    rel_ch = pat_ch{p(1)};

    patrates(p(2)) = sum(pat_ch96(p(1), rel_ch))/sum(N96(rel_ch));
    patcounts(p(2)) = sum(pat_ch96(p(1), rel_ch));
  end

  [patcounts_s, si] = sort(patcounts);
  patcounts_cs = cumsum(patcounts_s);

  spl_idx = find(patcounts_cs >= patcounts_cs(end)*0.5, 1);

  lo_pats = up(si(1:(spl_idx - 1)));
  hi_pats = up(si(spl_idx:end));

  fprintf('\n *** %s: %d/%d (%0.2f%%) *** \n\n', X.title{i}, spl_idx, length(up), 100*(1 - spl_idx/length(up)))

  % }}}

  %
  % save uniqued MAFs {{{

  %
  % hi {{{
  Ms_hi = Ms;
  Ms_hi.mut = reorder_struct(Ms_hi.mut, ismember(Ms_hi.mut.pat_idx, hi_pats));

  %per-channel exposures
  a_idx = ismember(a, hi_pats);
  Ms_hi.npat = sum(sparse(ch96(b(a_idx)), a(a_idx), 1), 2);

  [~, pui, puj] = unique([Ms_hi.mut.chr Ms_hi.mut.pos], 'rows');
  h = full(sparse(puj, nbidx(sub2ind(size(nbidx), Ms_hi.mut.ch_base, Ms_hi.mut.newbase_idx_rc)), 1, max(puj), 3));
  Ms_hi.mut = keep_fields(reorder_struct(Ms_hi.mut, pui), {'chr' 'pos' 'c1025' 'c65' 'c32' 'c512' 'is_coding' 'tier' 'gene_idx' 'ref_idx'});
  Ms_hi.mut.count_nb = h;
  Ms_hi.mut.count = sum(Ms_hi.mut.count_nb, 2);

  saveM(Ms_hi, ['mutation_data/MC3.align75.ICE_PoN-uniqued.mutratesplit_NMF/' X.names{i} '_hi.M'])

  fprintf('hi\n')
  count(Ms_hi.mut.count)

  % }}}

  %
  % lo {{{
  Ms_lo = Ms;
  Ms_lo.mut = reorder_struct(Ms_lo.mut, ismember(Ms_lo.mut.pat_idx, lo_pats));

  %per-channel exposures
  a_idx = ismember(a, lo_pats);
  Ms_lo.npat = sum(sparse(ch96(b(a_idx)), a(a_idx), 1), 2);

  [~, pui, puj] = unique([Ms_lo.mut.chr Ms_lo.mut.pos], 'rows');
  h = full(sparse(puj, nbidx(sub2ind(size(nbidx), Ms_lo.mut.ch_base, Ms_lo.mut.newbase_idx_rc)), 1, max(puj), 3));
  Ms_lo.mut = keep_fields(reorder_struct(Ms_lo.mut, pui), {'chr' 'pos' 'c1025' 'c65' 'c32' 'c512' 'is_coding' 'tier' 'gene_idx' 'ref_idx'});
  Ms_lo.mut.count_nb = h;
  Ms_lo.mut.count = sum(Ms_lo.mut.count_nb, 2);

  saveM(Ms_lo, ['mutation_data/MC3.align75.ICE_PoN-uniqued.mutratesplit_NMF/' X.names{i} '_lo.M'])

  fprintf('lo\n')
  count(Ms_lo.mut.count)

  % }}}

  % }}}
end

% }}}

%}}}

%}}}
