
% met_ind =  find_cell('adp_c', model.mets);
% met_ind =  find_cell('nadh_c', model.mets);
% e: Tolerance. It is like a lower level of the sum

load('varma_model.mat')

model = Varma_CharConcFluxMCAReadyFDP1;

%% Get all alternative Pools
e = 1e-6;

ForDaniel = {};
for i = 1:size(model.mets, 1)
    [coeffs_i, Pool_i] = alternatives(model, i, e);
    ForDaniel{i,1} = coeffs_i;
    ForDaniel{i,2} = Pool_i;
end

%% Get unique pools 
UniquePools = {};
k = 1;
for i = 1:size(model.mets, 1)
    if isempty(ForDaniel{i,2})
        continue        
    end
    contains_pool = 0;
    if ~isempty(UniquePools)
        for kp=1:size(UniquePools,1)
            if isequal(ForDaniel{i,2},UniquePools{kp,2})
                contains_pool = 1;
            end
        end
    end
    
    if contains_pool == 0 
        UniquePools{k,1} = ForDaniel{i,1};
        UniquePools{k,2} = ForDaniel{i,2};
        k = k + 1
    end
    
  
end

%% Create the conservation relation matrix 

L0_all = zeros(1, length(model.mets));
l = 1;
for i = 1:size(UniquePools, 1)
    % Unpack alternatives
    for k = 1:size(UniquePools{i,2},2)
        
        met_ids = UniquePools{i,2}{1,k};
        coeffs = UniquePools{i,1}(k);
        % TO MAKE THEM INTEGER 
        coeffs
        coeffs_int = coeffs/min(coeffs);
        coeffs_int
        met_idxes = [];
        for met_id = met_ids
            met_idx = find(ismember(model.mets, met_id));
            met_idxes = [met_idxes, met_idx];
        end
        
        L0_all(l,met_idxes) = coeffs_int;
        l = l + 1;
    end
end


%%
[ Abasis, L0i, lindeps]= getLinearIndependent(L0_all');

L0 = L0_all(L0i,:);

% Get depednent mets 
[L0_rref,depdenet_indexes] = rref(L0);
dep_met = model.mets(depdenet_indexes);

% Write a csv file 
L0_table = array2table(L0);
%L0_table.Properties.VariableNames = model.mets;
writetable(L0_table,'ecoli_nonlinear.csv')

% Print moieties 
for i = 1:size(L0,1)
    row = L0(i,:);
    met_ix = find(row);
    model.mets(met_ix)
end 
