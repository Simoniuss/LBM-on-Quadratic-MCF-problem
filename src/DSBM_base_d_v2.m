function [d, exitflag, v, theta, rho] = DSBM_base_d_v2(B_z, B_alpha, l, mu)

dSize = size(mu,1);

dVar = sdpvar(dSize,1);
vVar = sdpvar(1,1);

Constraints = [ (l >= vVar):'rho', (vVar >= B_z'*(dVar) - B_alpha'): 'theta'];
Objective = vVar + 0.1*0.5*norm(dVar)^2;

%options = sdpsettings('verbose', 0, 'solver','quadprog');
options = sdpsettings('verbose', 0, 'solver','mosek');
%options = sdpsettings('verbose', 0);

sol = optimize(Constraints, Objective, options);
d = value(dVar);
v = value(vVar);
theta = dual(Constraints('theta'));
rho = dual(Constraints('rho'));
exitflag = sol.problem;

end
