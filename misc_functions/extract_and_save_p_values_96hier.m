function extract_and_save_p_values_96hier(modeldir, sig)

F = [];
F.file = direc(['models_v3/' modeldir '/96_regr_hier/*/output/*.mat']);
F = parsein(F, 'file', '.*/(.*?)/output/(\d+)-([ACGT]\(.->.\)[ACGT]).*\.mat$', {'sig' 'ch96' 'context'});
F.ch96 = str2double(F.ch96);

if exist('sig', 'var'),
  F = reorder_struct(F, strcmp(F.sig, sig));
end

F = sort_struct(F, 'sig');

[~, sui] = unique(F.sig);

for x = [sui [sui(2:end) - 1; slength(F)]]',
  i = x(1); j = x(2);

  [L P] = process_run(reorder_struct(F, i:j));
  save(['models_v3/' modeldir '/96_regr_hier/' F.sig{i} '/loci_pvalues.mat'], 'L');
end
