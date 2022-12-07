%[Q, q, E, b, u] = loadMCF('MCF_1000.mat');
%[Q, q, E, b, u] = loadMCF('./data/netgen-20-1-0-b-b-ns.mat');
simpleProblem;

u(1) = 1;

xVar = sdpvar(size(q,1), 1);
Constraints = [ (E*xVar == b): 'equality', 0 <= xVar, xVar <= u ];
Objective = xVar' * diag(Q) * xVar + q'*xVar;

%options = sdpsettings('verbose', 0, 'solver','quadprog');
options = sdpsettings('verbose', 1, 'solver','mosek');
%options = sdpsettings('verbose', 1, 'solver','mosek', 'mosek.MSK_DPAR_INTPNT_QO_TOL_PFEAS',1e-12);
%options = sdpsettings('verbose', 0);

%savecplexlp(Constraints, Objective, 'probData');

sol = optimize(Constraints,Objective, options);
x = value(xVar);
fx = value(Objective);
mu = dual(Constraints('equality'));
exitflag = sol.problem;