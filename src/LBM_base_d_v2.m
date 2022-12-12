function [d, exitflag, theta] = LBM_base_d_v2(B_z, B_alpha, l, mu, fmu)

dSize = size(mu,1);

dVar = sdpvar(dSize,1);

Constraints = B_z'*(dVar) - B_alpha' <= l;
%Constraints = B_z'*(dVar) - B_alpha' + B_z'*mu - fmu  <= l;
Objective = 0.5*norm(dVar)^2;

%options = sdpsettings('verbose', 0, 'solver','quadprog');
options = sdpsettings('verbose', 0, 'solver','mosek');
%options = sdpsettings('verbose', 0);

sol = optimize(Constraints, Objective, options);
d = value(dVar);
theta = dual(Constraints);
exitflag = sol.problem;

end
