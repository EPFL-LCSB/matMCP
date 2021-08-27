function [coeffs, Pool] = alternatives(model, met_ind, e)
pools  = [];
coeffs = [];
k = 1;
[sol, a] = MCP(model, met_ind, e);
y = find(abs(sol.x(size(model.S, 1)+1:end)-1)<1e-7);
pools(1:length(y), k)  = y;
coeffs(1:length(y), k) = sol.x(y);
while ~isempty(sol.x)
    k = k+1;
    a.A(end+1, size(model.S, 1)+y) = 1;
    a.rhs(end+1)                 = length(y)-1;
    a.constraintType(end+1)      = {'<'};
    a.constraintNames(end+1)     = {'alt'};
    sol = solveTFAmodel(a);
    y = find(abs(sol.x(size(model.S, 1)+1:end)-1)<1e-7);
    pools(1:length(y), k)  = y;
    coeffs(1:length(y), k) = sol.x(y);
end
pools(:, end)  = [];
coeffs(:, end) = [];
Pool          = model.mets(pools);
end