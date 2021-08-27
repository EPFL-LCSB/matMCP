function [sol, a] = MCP(model, i, k)

Aeq = (model.S)';
Aeq = [Aeq, zeros(size(model.S, 2), size(model.S, 1))];
Aeq = [Aeq; [ones(1, size(model.S,1)), zeros(1, size(model.S, 1))]];
Aineq = [eye(size(model.S, 1)), -eye(size(model.S, 1))];
A     = [Aeq; Aineq];

lb = zeros(2*size(model.S,1),1);
lb(i) = k;
ub = ones(2*size(model.S,1),1);

f = zeros(2*size(model.S,1),1);
f(size(model.S,1)+1:end) = 1;

ctype = {};
ctype(1:size(model.S,1)) = {'C'};
ctype(size(model.S, 1)+1:2*size(model.S, 1)) = {'B'};
ctype = ctype';

rhs = [zeros(size(model.S, 2), 1); 1; zeros(size(model.S, 1), 1)];

colname = [model.mets; model.mets];
rowname = [model.rxns; {'a'}; model.mets];

a.A        = A;
a.rhs      = rhs;
a.varNames = colname;
a.vartypes = ctype;
a.var_lb   = lb;
a.var_ub   = ub;
a.f        = f;
a.objtype  = 1;
a.constraintType = {};

a.constraintType(1:size(model.S, 2)+1) = {'='};
a.constraintType(size(model.S, 2)+2:size(model.S, 2)+1+size(model.S, 1)) = {'<'};

a.constraintType = a.constraintType';
a.description    = model.description;

a.constraintNames = rowname;
sol               = solveTFAmodel(a);

end
