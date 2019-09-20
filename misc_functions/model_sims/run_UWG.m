function run_LNP(s)

if isdeployed,
  s = str2double(s);
end

load('ref/gene_list.align75_filtered.territories_UWG-Fg.mat', 'G')
N32 = sum(G.terr32);

M = loadM('mutation_data/MC3.align75.ICE_PoN-uniqued.M');
M.gene.gidx2 = listmap(M.gene.name, G.gene);
M.mut.gidx2 = M.gene.gidx2(M.mut.gene_idx);

%get full set of lambdas and their counts
lams = NaN(32, 3, slength(G));
for i = 1:32,
  idx = M.mut.c32 == i;

  terr = N32(i);

  for j = 1:3,
    ct = M.mut.count_nb(idx, j);
    lams(i, j, :) = sum(ct)/terr*G.Fg;
  end
end

terr_rep = permute(repmat(G.terr32, 1, 1, 3), [2 3 1]);

u = unique(G.terr32)';

%generate samples
full_samps = cell(100, 1);

pp = parpool(16);
parfor k = 1:100, 
  samps = cell(length(u), 1);
  for x = [u; 1:length(u)]
    i = x(1); j = x(2); 

    idx = terr_rep == i;

    samps{j} = poissrnd(repmat(lams(idx), i, 1));
  end
  full_samps{k} = cat(1, samps{:});
end
pp.delete

%concatenate and save
full_samps = sparse(cat(2, full_samps{:}));

save(sprintf('figures/model_sims/samps_UWG_%d.mat', s), 'full_samps', '-v7.3')
