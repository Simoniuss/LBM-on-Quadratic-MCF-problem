%[Q, q, E, b, u] = loadMCF('MCF_1000,mat');
simpleProblem;

xVar = sdpvar(size(q,1), 1);
Constraints = [ E*xVar == b, 0 <= xVar, xVar <= u ];
Objective = xVar' * diag(Q) * xVar + q'*xVar;

%options = sdpsettings('verbose', 0, 'solver','quadprog');
options = sdpsettings('verbose', 1, 'solver','mosek');
%options = sdpsettings('verbose', 0);

sol = optimize(Constraints,Objective, options);
x = value(xVar);
fx = value(Objective);
exitflag = sol.problem;