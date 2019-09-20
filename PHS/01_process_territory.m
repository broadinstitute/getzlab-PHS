% NOTE: the raw covariate data processed by this script into a format 
%       parseable by the regression models is quite large (>100 GB); as such, we are only
%       including the processed versions of these covariates (i.e., the output
%       of this script).  Raw covariate data available on request.

%
% 0. mask areas of low mappability in target list {{{
clear

%load target list
T = load_struct('ref/gencode/target_list.v1.txt');
T = makeapn(T);

%
%get mappability across targets
tmap = cell(slength(T), 1);

parfor i = 1:length(tmap),
  fwb = org.broadinstitute.cga.tools.seq.FixedWidthBinary('ref/map/align75.fwb');
  tmap{i} = fwb.get(T.chr(i), T.start(i), T.end(i));
  fwb.close();
end

T.tmap = cellfun(@(x) 10.^(5.64/256*double(x)), tmap, 'unif', 0);

%let's be very conservative: restrict ourselves to completely unique 75mers
T.mask = cellfun(@(x) x == 1, T.tmap, 'unif', 0);

%
%segment targets according to mask
T.ustarts = cell(slength(T), 1);
T.uends = cell(slength(T), 1);
for i = 1:slength(T),
  d = diff([0; T.mask{i}; 0]);
  T.ustarts{i} = find(d == 1) + T.start(i) - 1;
  T.uends{i} = find(d == -1) + T.start(i) - 2;

  %conservatively pad segments, such that no gap can be < 24 bases long)
  T.ustarts{i}(2:end) = T.ustarts{i}(2:end) + 12;
  T.uends{i}(1:(end - 1)) = T.uends{i}(1:(end - 1)) - 12;
end

%
%generate new target list according to these unique segments
ntargs = sum(cellfun(@length, T.ustarts));

Tu = [];
Tu.start = NaN(ntargs, 1);
Tu.end = NaN(ntargs, 1);

Tu.chr = NaN(ntargs, 1);
Tu.gene = cell(ntargs, 1);

c = 1;
for i = 1:slength(T),
  l = length(T.ustarts{i});
  if l == 0, continue; end

  r = c:(c + l - 1);

  Tu.start(r) = T.ustarts{i};
  Tu.end(r) = T.uends{i};

  Tu.chr(r) = T.chr(i);
  Tu.gene(r) = T.gene(i);

  c = c + l;
end

%
%save
Tu = reorder_struct(Tu, Tu.end > Tu.start);

save_struct(order_fields(Tu, {'gene' 'chr' 'start' 'end'}), 'ref/target_list.align75_filtered.txt')

% }}}

%
% 1. generate covariate directories % {{{
clear

% copy context binaries to shared memory for quick access
% this needs at least 3 GB of free memory.
% cp ref/pentamer/chr*.mat /dev/shm

context_dir = '/dev/shm/';

%make for all coding
T = load_struct('ref/target_list.align75_filtered.txt');
T = makeapn(T);

pp = parpool(12);

%
%coarse covariates {{{
C = [];
C.file = direc('ref/raw_transformed_covars_noZ/*.mat');
C.covar = regexprep(C.file, '.*/(.*)\.mat$', '$1');
C.outdir = strcat('ref/cod_mutsig_transformed_covars.align75_filtered/', C.covar);

parfor i = 1:slength(C),
  covar_index(context_dir, T, C.file{i}, C.outdir{i}, 1);
end

% }}}

%
%XR-seq {{{
F = [];
F.file = direc('ref/XRseq/*.fwb');
%NHF1, XPC, CSBC
F = parsein(F, 'file', '.*7941_(.*)(64|CPD).*_Merged_(MINUS|PLUS).*', {'cellline' 'lesion' 'strand'});
F.cellline(ismember(F.cellline, {'' 'PP'})) = {'NHF1'};

F.outdir = strcat('ref/cod_XR_covars.align75_filtered/', F.cellline, '_', F.lesion, '_', F.strand);

parfor i = 1:slength(F),
  covar_index(context_dir, T, F.file{i}, F.outdir{i}, 1);
end

% }}}

%
%nucleosome positions {{{

%see db/hg19/nuc/ for FWB conversions
F = [];
F.file = direc('ref/nuc/*.fwb');
F.cellline = regexprep(F.file, '.*wgEncodeSydhNsome(.*)Sig.*', '$1');
F.outdir = strcat('ref/cod_nuc_covars.align75_filtered/', F.cellline);

parfor i = 1:slength(F),
  covar_index(context_dir, T, F.file{i}, F.outdir{i}, 1);
end

% }}}

%
%DNAse HS sites {{{

%see db/hg19/dnasehs/ for bedfile conversions
covar_index(context_dir, T, 'ref/dnasehs/wgEncodeRegDnaseClusteredV3.scoreonly.zerofilled.mat', 'ref/cod_dnasehs_covars.zerofilled.align75_filtered', 1);

% mkdir ref/cod_dnasehs_covars.zerofilled.align75_filtered/v3
% mv ref/cod_dnasehs_covars.zerofilled.align75_filtered/*.mat ref/cod_dnasehs_covars.zerofilled.align75_filtered/v3

% }}}

%
%pentamers (16) {{{

%
%we first have to make an FWB of strand-collapsed contexts
% (i.e. c512)
%do this in ~/jhcga/db/hg19/context512

%
%now index WRT this track

covar_index(context_dir, T, 'ref/context512/all.fwb', 'ref/cod_c512_covars.align75_filtered');

%}}}

%
% APOBEC hot hairpins {{{

load('ref/apobec/all_TpC_sites.v2.0.mat', 'S');

is_coding = map_mutations_to_targets_fast([double(S.chr) double(S.pos)], [T.chr T.start T.end]) > 0;

S = reorder_struct(S, is_coding);
S.chr = double(S.chr);
S.start = double(S.pos);
S.end = S.start;
S.covar = log10(double(S.relrate_exp));

Q = keep_fields(S, {'chr' 'start' 'end' 'covar'});

save('ref/raw_covars/APOBEC.lawrence.hairpins_v2.mat', 'Q')

covar_index(context_dir, T, 'ref/raw_covars/APOBEC.lawrence.hairpins_v2.mat', 'ref/cod_APOBEC.lawrence.hairpin_covars.align75_filtered/v2');

% }}}

%cleanup
pp.delete
% rm -f /dev/shm/chr*.mat

% }}}

%
% 2. generate territory table {{{
clear

% this needs at least 3 GB of free memory.
% cp ref/pentamer/chr*.mat /dev/shm

T = load_struct('ref/target_list.align75_filtered.txt');
T = makeapn(T);
T = sort_struct(T, {'chr' 'start' 'end'});

[cu, cui] = unique(T.chr);
[~, gui, T.gidx] = unique(T.gene);

h = zeros(1024, 1);
for x = [cui [cui(2:end) - 1; slength(T)] (1:24)']',
  i = x(1); j = x(2); k = x(3);

  load(sprintf('/dev/shm/chr%d.mat', k))

  is_coding = false(size(categ)); 
  for t = i:j,
    is_coding(T.start(t):T.end(t)) = true;
  end

  h = h + histc(categ(is_coding), 1:1024);
end

save('ref/pentamer_territories.align75_filtered.mat', 'h');

%
% 2.1 also compute generate gene-specific territories {{{

G = [];
G.gene = T.gene(gui);
G.terr = NaN(slength(G), 1024);

for x = [cui [cui(2:end) - 1; slength(T)] cu]',
  i = x(1); j = x(2); k = x(3);

  load(sprintf('/dev/shm/chr%d.mat', k))

  guii = gui(T.chr(gui) == k);

  is_gene = false(size(categ));
  for y = [guii [guii(2:end) - 1; j]]',
    u = y(1); v = y(2);

    for t = u:v,
      is_gene(T.start(t):T.end(t)) = true;
    end

    chridx = find(is_gene(T.start(u):T.end(v))) + T.start(u) - 1;
    G.terr(T.gidx(u), :) = histc(categ(chridx), 1:1024);

    for t = u:v,
      is_gene(T.start(t):T.end(t)) = false;
    end
  end

  fprintf('%d ', k)
end

%collapse to trimers
c1024 = dec2base(0:1023, 4) - 48;
c64 = c1024(:, [1 3 2])*[16 4 1]' + 1;

G.terr64 = NaN(slength(G), 64);

cols = accumarray(c64, 1:1024, [], @(x) {x});
for i = 1:length(cols),
  G.terr64(:, i) = sum(G.terr(:, cols{i}), 2);
end

%strand collapse
c64 = dec2base(0:63, 4) - 48;
c64 = [(1:64)' sum(bsxfun(@times, 3 - c64(:, [1 3 2]), [4^2 4 1]), 2) + 1]; 
c32map(c64(33:end, 2)) = c64(33:end, 1);

G.terr32 = G.terr64(:, 1:32) + G.terr64(:, c32map);

%save
save('ref/gene_list.align75_filtered.territories.mat', 'G')

% }}}

% rm -f /dev/shm/chr*.mat

%}}}
