%[Q, q, E, b, u] = loadMCF('MCF_1000.mat');
simpleProblem;

u(1) = 1;

xVar = sdpvar(size(q,1), 1);
Constraints = [ E*xVar == b, 0 <= xVar, xVar <= u ];
Objective = xVar' * diag(Q) * xVar + q'*xVar;

%options = sdpsettings('verbose', 0, 'solver','quadprog');
options = sdpsettings('verbose', 1, 'solver','mosek');
%options = sdpsettings('verbose', 1, 'solver','mosek', 'mosek.MSK_DPAR_INTPNT_QO_TOL_PFEAS',1e-12);
%options = sdpsettings('verbose', 0);

%savecplexlp(Constraints, Objective, 'probData');

sol = optimize(Constraints,Objective, options);
x = value(xVar);
fx = value(Objective);
exitflag = sol.problem;