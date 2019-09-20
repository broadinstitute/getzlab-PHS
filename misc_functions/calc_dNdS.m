function G = calc_dNdS(M, Mi, C, effect_terrs_by_gene)

%lookup tables for ch96 -> c32 x 3
tmp = dec2base(0:95, 4) - 48;
c96_32map = [tmp(:, 1:2)*[4 1]' > 2 tmp(:, 3:4)]*[16 4 1]' + 1;
c96_bc = mod(ceil((1:96)/16) - 1, 3)' + 1;

%total number of covered trimers in exome
N32 = sum(sum(effect_terrs_by_gene(:, :, 1, :), 2), 4);
N96 = N32(c96_32map);

%rate of each base change
r96 = histc(M.ch96, 1:96)./N96;
r96_ch = full(sparse(c96_32map, c96_bc, r96));

%rate of indels
r_ind = slength(Mi)/sum(N32);

G = C.gene;
G.tier = get_tiers(G.name, 'ref/gencode/gene_tiers.v10.3.mat');
G.tier(isnan(G.tier)) = 4;

%arrays for storing d[E]/dS

%all nonsyn:syn
dNdS_95CI = cell(slength(G), 1);
dN_95CI = cell(slength(G), 1);
dS_95CI = cell(slength(G), 1);

%truncating:syn
dTdS_95CI = cell(slength(G), 1);
dT_95CI = cell(slength(G), 1);

f_gt_05 = NaN(slength(G), 1);
f_in_neu = NaN(slength(G), 1);

parfor g = 1:size(effect_terrs_by_gene, 4)
  %expected number of mutations of each effect in each channel (Poisson means)
  N_syn = bsxfun(@times, r96_ch, squeeze(effect_terrs_by_gene(:, 1, :, g)));
  N_non = bsxfun(@times, r96_ch, squeeze(sum(effect_terrs_by_gene(:, 2:4, :, g), 2)));
  N_tru = bsxfun(@times, r96_ch, squeeze(sum(effect_terrs_by_gene(:, 3:4, :, g), 2)));

  %overall observed number of mutations of each effect in each channel
  idx = M.gene_idx_lm == g & M.effect_idx == 5;
  n_syn = full(sparse(M.c32(idx), M.nb_rc(idx), 1));

  idx = M.gene_idx_lm == g & ismember(M.effect_idx, [1 3 4]);
  n_non = full(sparse(M.c32(idx), M.nb_rc(idx), 1));

  idx = M.gene_idx_lm == g & ismember(M.effect_idx, [3 4]);
  n_tru = full(sparse(M.c32(idx), M.nb_rc(idx), 1));

  %totals of expected/observed (for input into likelihood)
  ss = fullsum(n_syn);
  sS = fullsum(N_syn);
  sn = fullsum(n_non) + nnz(Mi.gene_idx_lm == g);
  sN = fullsum(N_non) + r_ind*fullsum(effect_terrs_by_gene(:, :, 1, g));
  st = fullsum(n_tru) + nnz(Mi.gene_idx_lm == g);
  sT = fullsum(N_tru) + r_ind*fullsum(effect_terrs_by_gene(:, :, 1, g));

  %dN/dS is quotient of gamma random variables (FYI, this is the beta prime distribution)
  gr_n = gamrnd(sn + 1, 1/sN, 1e4, 1);
  gr_s = gamrnd(ss + 1, 1/sS, 1e4, 1);

  dNdS = gr_n./gr_s;
  dNdS_95CI{g} = [quantile(dNdS, [0.025 0.975]) mean(dNdS)];

  f_gt_05(g) = mean(log2(dNdS) > 0.5);
  f_in_neu(g) = mean(log2(dNdS) > -0.23 & log2(dNdS) < 0.13);

  dN_95CI{g} = [quantile(gr_n, [0.025 0.975]) mean(gr_n)];
  dS_95CI{g} = [quantile(gr_s, [0.025 0.975]) mean(gr_s)];

  %dT/dS
  gr_t = gamrnd(st + 1, 1/sT, 1e4, 1);
  dT_95CI{g} = [quantile(gr_t, [0.025 0.975]) mean(gr_t)];
  dTdS_95CI{g} = [quantile(gr_t./gr_s, [0.025 0.975]) mean(gr_t./gr_s)];
end

G.dNdS_95CI = cat(1, dNdS_95CI{:});
G.dN_95CI = cat(1, dN_95CI{:});
G.dS_95CI = cat(1, dS_95CI{:});

G.dTdS_95CI = cat(1, dTdS_95CI{:});
G.dT_95CI = cat(1, dT_95CI{:});

G.f_gt_05 = f_gt_05;
G.f_in_neu = f_in_neu;
