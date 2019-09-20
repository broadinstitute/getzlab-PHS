function X = generate_signature_subcohort_definitions(mafpath)

X = [];
X.mafpath = direc(mafpath);
X.names = regexprep(X.mafpath, '.*/([A-Z_]+)\.M$', '$1');
X.title = mapacross(X.names, ...
{'APOBEC' ...
'ESO' ...
'LUNG' ...
'MSI' ...
'POLE' ...
'POLE_MSI' ...
'UV' ...
'VANILLA'}, ...
{'APOBEC' ...
'Eso.' ...
'Smoking' ...
'MSI' ...
'POLE' ...
'POLE+MSI' ...
'UV' ...
'CpG (Aging)'});

X.context_regex = {{'C in .T_[ACT].'}; ... %APO
                   {'A in .A_..' 'A in .._G.'}; ... %ESO
                   {'C in .._..'}; ... %LUNG
                   {'C in .C_..' 'C in .G_..'}; ... %MSI
                   {'C in .T_T.' 'C in .T_[GT].' 'A in .A_A.'}; ... %POLE
                   {'C in .._T.' 'C in .T_[GT].'}; ... %POLE_MSI
                   {'C in .[TC]_..'}; ... %UV
                   {'C in .._G.'}}; %VANILLA
X.newbase = {{[1 2 3]}; ... %APO 
             {1 2}; ... %ESO 
             {1}; ... %LUNG
             {1 3}; ... %MSI
             {1 3 1}; ... %POLE
             {1 3}; ... %POLE_MSI
             {3}; ... %UV
             {3}}; ... %VANILLA
X.npat = NaN(slength(X), 1);

%define channels of interest
C = load_struct('ref/context_1025_categs.txt');
C = makeapn(C);

X.output = cell(slength(X), 1);

for i = 1:slength(X),
  X.output{i} = cell(size(X.context_regex{i}));
  for j = 1:length(X.context_regex{i})
    c512s = find(grepm(X.context_regex{i}{j}, C.name));
    [c nb] = meshgrid(c512s, X.newbase{i}{j});

    X.output{i}{j} = [c(:) nb(:)];
  end

  M = loadM(X.mafpath{i});
  X.npat(i) = full(max(M.npat)); %if we are reading in results produced from NMF, number of patients will vary based on channel ... for summary statistic here, just take the max.
end
