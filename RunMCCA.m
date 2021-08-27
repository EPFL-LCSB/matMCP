%restoredefaultpath

% load('/Users/georgiosfengos/Dropbox/SharedFolders/EPFL-LCSB/People/SharedDanielW/moieties/inputs.mat')
% load('varma_model.mat')

%addpath(genpath('/Users/georgiosfengos/GIT_Folders/FBA_Toolboxes'))
%addpath(genpath('/Users/georgiosfengos/Applications/IBM/ILOG/CPLEX_Studio1271'))

% met_ind =  find_cell('adp_c', model.mets);
% met_ind =  find_cell('nadh_c', model.mets);
% e: Tolerance. It is like a lower level of the sum
e = 1e-6;
ForDaniel = {};
for i = 1:size(model.mets, 1)
    [coeffs_i, Pool_i] = alternatives(model, i, e);
    ForDaniel{i,1} = coeffs_i;
    ForDaniel{i,2} = Pool_i;
end