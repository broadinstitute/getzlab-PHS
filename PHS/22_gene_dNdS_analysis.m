%for this, we will use the full MAF/territory
%we will also include indels.

%
%load target list/coverage models {{{

load('ref/gencode/coverage_models.v3.mat')

%target list
T = load_struct('ref/gencode/target_list.v1.txt');
T = makeapn(T);

%splice sites are +/- 2 from exon boundaries -- we need to include these
T.start = T.start - 2;
T.end = T.end + 2;

T = sort_struct(T, {'chr' 'start' 'end'});

%}}}

%
%load all sSNVs {{{
M = loadM('mutation_data/MC3.M');

%
%map mutations to MutSig geneset
M.mut = sort_struct(M.mut, {'chr' 'pos'});

%add back intronic splice variants
targ = map_mutations_to_targets_fast([M.mut.chr M.mut.pos], [T.chr T.start T.end]);
T.gene_idx = listmap(T.gene, M.gene.name);

splidx = targ > 0 & ~M.mut.is_coding & M.mut.effect_idx == 4;

M.mut.gene_idx(splidx) = T.gene_idx(targ(splidx));

%map to MutSig genes
M.mut.gene_idx_lm = listmap(M.gene.name(M.mut.gene_idx), C.gene.name); M.mut = reorder_struct(M.mut, ~isnan(M.mut.gene_idx_lm));

%collapse newbase strand
M.mut.newbase_idx_rc = M.mut.newbase_idx;
idx = M.mut.c1025 > 512;
M.mut.newbase_idx_rc(idx) = 5 - M.mut.newbase_idx(idx); 

ch_base = dec2base(M.mut.c32 - 1, 4, 3) - 47; M.mut.ch_base = ch_base(:, 1);
nbidx = [0 1 2 3; ...
         1 0 2 3; ...
         1 2 0 3; ...
         1 2 3 0];
M.mut.nb_rc = nbidx(sub2ind(size(nbidx), M.mut.ch_base, M.mut.newbase_idx_rc));

%add binomial q-values (to remove "significant" hotspots to see how dN/dS changes)
load('bino/loci_pvalues.mat', 'L')

M.mut = multimapinto(M.mut, L, {'chr' 'pos' 'ch1536'}, {'chr' 'pos' 'ch1536'}, {'p_bino' 'q_bino'});
M.mut.sig_bino = M.mut.q_bino <= 0.25;

%}}}

%
%load indels {{{
Mi = loadM('mutation_data/MC3_indels.M');

%map mutations to MutSig geneset
Mi.mut.gene_idx_lm = listmap(Mi.gene.name(Mi.mut.gene_idx), C.gene.name); Mi.mut = reorder_struct(Mi.mut, ~isnan(Mi.mut.gene_idx_lm));

%}}}

%
%calculate overall rates {{{

%load territories
load('ref/trimer_coding_effect_territories_by_gene.mat')

%lookup tables for ch96 -> c32 x 3
tmp = dec2base(0:95, 4) - 48;
c96_32map = [tmp(:, 1:2)*[4 1]' > 2 tmp(:, 3:4)]*[16 4 1]' + 1;
c96_bc = mod(ceil((1:96)/16) - 1, 3)' + 1;

%total number of covered trimers in exome
N32 = sum(sum(effect_terrs_by_gene(:, :, 1, :), 2), 4);
N96 = N32(c96_32map);

%rate of each base change
r96 = histc(M.mut.ch96, 1:96)./N96;
r96_ch = full(sparse(c96_32map, c96_bc, r96));

%rate of indels
r_ind = slength(Mi.mut)/sum(N32);

%}}}

%
%compute dN/dS for each gene, for each ttype {{{
G = C.gene;
G.tier = get_tiers(G.name, 'gencode/gene_tiers.v10.3.mat');
G.tier(isnan(G.tier)) = 4;

%arrays for storing d[E]/dS

%all nonsyn:syn
dNdS_95CI = cell(slength(G), 1);
dN_95CI = cell(slength(G), 1);
dS_95CI = cell(slength(G), 1);
f_gt_05 = NaN(slength(G), 1);
f_in_neu = NaN(slength(G), 1);
dNdS_95CI_nohot = cell(slength(G), 1);

%missense:syn
dMdS_95CI = cell(slength(G), 1);
dM_95CI = cell(slength(G), 1);
dMdS_95CI_nohot = cell(slength(G), 1);

%truncating:syn
dTdS_95CI = cell(slength(G), 1);
dT_95CI = cell(slength(G), 1);
dTdS_95CI_nohot = cell(slength(G), 1);

%for dN/dS, we also compute histograms
dNdS_hist = cell(slength(G), 1);
dNdS_hist_fine = cell(slength(G), 1);

hist_range = linspace(0, 35, 1e3);
hist_range_fine = linspace(0, 2, 1e3);

pp = parpool(20);
parfor g = 1:size(effect_terrs_by_gene, 4)
  %expected number of mutations of each effect in each channel (Poisson means)
  N_syn = bsxfun(@times, r96_ch, squeeze(effect_terrs_by_gene(:, 1, :, g)));
  N_non = bsxfun(@times, r96_ch, squeeze(sum(effect_terrs_by_gene(:, 2:4, :, g), 2)));
  N_mis = bsxfun(@times, r96_ch, squeeze(effect_terrs_by_gene(:, 2, :, g)));
  N_tru = bsxfun(@times, r96_ch, squeeze(sum(effect_terrs_by_gene(:, 3:4, :, g), 2)));

  %overall observed number of mutations of each effect in each channel
  idx = M.mut.gene_idx_lm == g & M.mut.effect_idx == 5;
  n_syn = full(sparse(M.mut.c32(idx), M.mut.nb_rc(idx), 1));

  idx = M.mut.gene_idx_lm == g & ismember(M.mut.effect_idx, [1 3 4]);
  n_non = full(sparse(M.mut.c32(idx), M.mut.nb_rc(idx), 1));

  misidx = M.mut.gene_idx_lm == g & M.mut.effect_idx == 1;
  n_mis = full(sparse(M.mut.c32(misidx), M.mut.nb_rc(misidx), 1));

  idx = M.mut.gene_idx_lm == g & ismember(M.mut.effect_idx, [3 4]);
  n_tru = full(sparse(M.mut.c32(idx), M.mut.nb_rc(idx), 1));

  %totals of expected/observed (for input into likelihood)
  ss = fullsum(n_syn);
  sS = fullsum(N_syn);
  sn = fullsum(n_non) + nnz(Mi.mut.gene_idx_lm == g);
  sN = fullsum(N_non) + r_ind*fullsum(effect_terrs_by_gene(:, :, 1, g));
  sm = fullsum(n_mis);
  sM = fullsum(N_mis);
  st = fullsum(n_tru) + nnz(Mi.mut.gene_idx_lm == g);
  sT = fullsum(N_tru) + r_ind*fullsum(effect_terrs_by_gene(:, :, 1, g));

  %dN/dS is quotient of gamma random variables
  gr_n = gamrnd(sn + 1, 1/sN, 1e4, 1);
  gr_s = gamrnd(ss + 1, 1/sS, 1e4, 1);

  dNdS = gr_n./gr_s;
  dNdS_95CI{g} = [quantile(dNdS, [0.025 0.975]) mean(dNdS)];

  f_gt_05(g) = mean(log2(dNdS) > 0.5);
  f_in_neu(g) = mean(log2(dNdS) > -0.23 & log2(dNdS) < 0.13);

  dN_95CI{g} = [quantile(gr_n, [0.025 0.975]) mean(gr_n)];
  dS_95CI{g} = [quantile(gr_s, [0.025 0.975]) mean(gr_s)];

  dNdS_hist{g} = histc(dNdS, hist_range);
  dNdS_hist_fine{g} = [histc(dNdS, hist_range_fine); nnz(dNdS > 2)];

  %dT/dS
  gr_t = gamrnd(st + 1, 1/sT, 1e4, 1);
  dT_95CI{g} = [quantile(gr_t, [0.025 0.975]) mean(gr_t)];
  dTdS_95CI{g} = [quantile(gr_t./gr_s, [0.025 0.975]) mean(gr_t./gr_s)];

  %for genes that contain missense hotspots (very liberally defined)
  %we calculate the drop in dM/dS after removing the hotspots

  if any(M.mut.sig_bino(misidx)),
    idx = M.mut.gene_idx_lm == g & M.mut.effect_idx == 1 & ~M.mut.sig_bino;
    n_mis = full(sparse(M.mut.c32(idx), M.mut.nb_rc(idx), 1));

    sm_nohot = fullsum(n_mis);

    %dM/dS excluding hotspots
    gr = gamrnd(sm_nohot + 1, 1/sM, 1e4, 1)./gamrnd(ss + 1, 1/sS, 1e4, 1); 
    dMdS_95CI_nohot{g} = [quantile(gr, [0.025 0.975]) mean(gr)];

    %we only care about the following when evaluating potential TP candidates.

    %dM/dS including hotspots
    gr = gamrnd(sm + 1, 1/sM, 1e4, 1)./gamrnd(ss + 1, 1/sS, 1e4, 1); 
    dMdS_95CI{g} = [quantile(gr, [0.025 0.975]) mean(gr)];
  end
end

G.dNdS_95CI = cat(1, dNdS_95CI{:});
G.dN_95CI = cat(1, dN_95CI{:});
G.dS_95CI = cat(1, dS_95CI{:});
G.dT_95CI = cat(1, dT_95CI{:});
G.dNdS_hist = cat(2, dNdS_hist{:})';
G.dNdS_hist_fine = cat(2, dNdS_hist_fine{:})';

G.dTdS_95CI = cat(1, dTdS_95CI{:});

save('G_dNdS.mat', 'G')

%true positives
idx = ~cellfun(@isempty, dMdS_95CI_nohot);
G_tp = reorder_struct(G, idx);
G_tp.dMdS_95CI = cat(1, dMdS_95CI{:});
G_tp.dMdS_95CI_nohot = cat(1, dMdS_95CI_nohot{:});
G_tp.dMdS_change = log2(G_tp.dMdS_95CI(:, 3)) - log2(G_tp.dMdS_95CI_nohot(:, 3));

save('G_dNdS_truepositives.mat', 'G_tp')

%
% per-ttypes {{{

%above code broken out into calc_dNdS.m
Gc = cell(slength(M.ttype), 1);


for i = 1:slength(M.ttype),
  Mt = reorder_struct(M.mut, M.mut.ttype_idx == i);
  Mit = reorder_struct(Mi.mut, Mi.mut.ttype_idx == i);

  Gc{i} = calc_dNdS(Mt, Mit, C, effect_terrs_by_gene);
end

Gtt = Gc{1};
Gtt.dNdS_95CI_ttype = NaN(slength(Gtt), 3, slength(M.ttype));
Gtt.dN_95CI_ttype = NaN(slength(Gtt), 3, slength(M.ttype));
Gtt.dS_95CI_ttype = NaN(slength(Gtt), 3, slength(M.ttype));
Gtt.dTdS_95CI_ttype = NaN(slength(Gtt), 3, slength(M.ttype));
Gtt.dT_95CI_ttype = NaN(slength(Gtt), 3, slength(M.ttype));
Gtt.f_gt_05_ttype = NaN(slength(Gtt), slength(M.ttype));
Gtt.f_in_neu_ttype = NaN(slength(Gtt), slength(M.ttype));

for i = 1:slength(M.ttype),
  Gtt.dNdS_95CI_ttype(:, :, i) = Gc{i}.dNdS_95CI;
  Gtt.dN_95CI_ttype(:, :, i) = Gc{i}.dN_95CI;
  Gtt.dS_95CI_ttype(:, :, i) = Gc{i}.dS_95CI;
  Gtt.dTdS_95CI_ttype(:, :, i) = Gc{i}.dTdS_95CI;
  Gtt.dT_95CI_ttype(:, :, i) = Gc{i}.dT_95CI;
  Gtt.f_gt_05_ttype(:, i) = Gc{i}.f_gt_05;
  Gtt.f_in_neu_ttype(:, i) = Gc{i}.f_in_neu;
end

save('G_dNdS_perttype.mat', 'Gtt')

%}}}

%}}}

%
% Table S3 {{{
clear

load('G_dNdS.mat', 'G')

idx1 = G.dTdS_95CI(:, 1) >= 2.5;
dif = diff(log2(G.dT_95CI(:, 1:2)), [], 2);

Gts = reorder_struct(G, idx1 & dif < 0.95);

% map in CGC status
load('ref/CGC/v1.1.mat', 'C')
Gts.is_KCG = nansub(C.is_KCG, listmap(Gts.name, C.GeneSymbol), false);

% save table
Gts.CI95_lo = cellsprintf('%0.2f', Gts.dTdS_95CI(:, 1));
Gts.CI95_lo_num = Gts.dTdS_95CI(:, 1);
save_struct(keep_fields(sort_struct(Gts, {'is_KCG' 'tier' 'CI95_lo_num'}, [-1 1 1]), {'name' 'tier' 'is_KCG' 'CI95_lo'}), 'tables/TS_genes_supp.tsv')

%}}}

%
% Table S1 {{{

clear

load('G_dNdS.mat', 'G')
load('G_dNdS_truepositives.mat', 'G_tp')
load('G_dNdS_perttype.mat', 'Gtt');

% index genes under positive selection
G.dMdS_change = mapacross(G.name, G_tp.name, G_tp.dMdS_change);

G.pos_idx = ismember(G.name, G_tp.name(G_tp.tier < 4 & G_tp.dMdS_change > log2(1.05) & G_tp.dTdS_95CI(:, 1) < 1.5));

% index genes under neutral selection (Prob[0.8 < dN/dS < 1.2] > 0.95)
hist_range_fine = linspace(0, 2, 1e3);

rng = [find(hist_range_fine > 0.8, 1) find(hist_range_fine > 1.2, 1)];
G.dNdS_p = sum(G.dNdS_hist_fine(:, rng(1):rng(2)), 2)/10000;
G.neg_idx = G.dNdS_p > 0.95;

%map in q-values from methods
M = loadM('mutation_data/MC3.align75.ICE_PoN-pvalues_a212692_v2.M');

%minimum q-value for all sites in each gene
Gtmp = [];
Gtmp.q = accumarray(M.mut.gene_idx, M.mut.q(:, 3), [], @min);
Gtmp.q_unif = accumarray(M.mut.gene_idx, M.mut.q_unif, [], @min);
Gtmp.q_uwg = accumarray(M.mut.gene_idx, M.mut.q_uwg, [], @min);
Gtmp.q_nb = accumarray(M.mut.gene_idx, M.mut.q_nb, [], @min);
Gtmp.gene = M.gene.name(1:slength(Gtmp));

G = mapinto(G, Gtmp, 'name', 'gene', {'q' 'q_unif' 'q_uwg' 'q_nb'});

G.sig_UWG = G.q_uwg < 0.1;
G.sig_UP = G.q_unif < 0.1;
G.sig_NB = G.q_nb < 0.1;
G.sig_LNP = G.q < 0.1;

%
% save FP table (sheet 1) {{{

G_fp = reorder_struct(G, G.neg_idx);
G_fp = sort_struct(G_fp, 'dNdS_p', -1);

G_fp.dNdS_mean = G_fp.dNdS_95CI(:, 3);
G_fp.dNdS_95CI = strcat('[', cellsprintf('%0.2f', G_fp.dNdS_95CI(:, 1)), {' '}, cellsprintf('%0.2f', G_fp.dNdS_95CI(:, 2)), ']');

save_struct(rename_field(keep_fields(G_fp, {'name' 'dNdS_p' 'dNdS_mean' 'dNdS_95CI' 'sig_UWG' 'sig_UP' 'sig_NB' 'sig_LNP'}), 'dNdS_p', 'neu_prob'), 'tables/FP_truth_set.tsv')

% }}}

%
% save TP table (sheet 2) {{{

G_tp = reorder_struct(G, G.pos_idx);
G_tp = sort_struct(G_tp, 'dMdS_change', -1);

G_tp.dNdS_mean = G_tp.dNdS_95CI(:, 3);
G_tp.dNdS_95CI = strcat('[', cellsprintf('%0.2f', G_tp.dNdS_95CI(:, 1)), {' '}, cellsprintf('%0.2f', G_tp.dNdS_95CI(:, 2)), ']');

save_struct(rename_field(keep_fields(G_tp, {'name' 'dMdS_change' 'dNdS_mean' 'dNdS_95CI' 'sig_UWG' 'sig_UP' 'sig_NB' 'sig_LNP'}), 'dMdS_change', 'dNdS_change'), 'tables/TP_truth_set.tsv')

% }}}

%
% save handpicked TP table (sheet 3) {{{

tp_list = {'AKT1' ...
'B2M' ...
'BCOR' ...
'BRAF' ...
'CDKN2A' ...
'CHD4' ...
'CIC' ...
'CTNNB1' ...
'EGFR' 'EP300' ...
'ERBB2' ...
'ERBB3' ...
'FBXW7' ...
'FGFR2' ...
'FGFR3' ...
'GNA11' ...
'GNAQ' ...
'GNAS' ...
'HRAS' ...
'IDH1' ...
'IDH2' ...
'KEAP1' ...
'KIT' ...
'KRAS' ...
'MAP2K1' ...
'MAPK1' ...
'MAX' 'MED12' ...
'MET' ...
'MTOR' ...
'NFE2L2' ...
'NOTCH1' ...
'NRAS' ...
'PIK3CA' 'PIK3CB' ...
'PIK3R1' ...
'POLE' ...
'PPP2R1A' ...
'RAC1' ...
'SF3B1' ...
'SMAD4' ...
'SMARCA4' ...
'SPOP' ...
'TP53' ...
'U2AF1' ...
'ACTB' ...
'CDK4' ...
'ERCC2' ...
'GNA13' ...
'KLF5' ...
'MYCN' ...
'SOX17' ...
'STK19'};

G_tp = reorder_struct(G, ismember(G.name, tp_list));

G_tp.dNdS_mean = G_tp.dNdS_95CI(:, 3);
G_tp.dNdS_95CI = strcat('[', cellsprintf('%0.2f', G_tp.dNdS_95CI(:, 1)), {' '}, cellsprintf('%0.2f', G_tp.dNdS_95CI(:, 2)), ']');

save_struct(keep_fields(G_tp, {'name' 'dNdS_mean' 'dNdS_95CI' 'sig_UWG' 'sig_UP' 'sig_NB' 'sig_LNP'}), 'tables/TP_handpicked_truth_set.tsv')

% }}}

% }}}
