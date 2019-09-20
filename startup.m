addpath('./funcs')

[~, p] = unix('find funcs misc_functions LNP -type d -not -path ''*/\.*''');
p = strsplit(p, '\n');

for i = p, addpath(decell(i)); end
