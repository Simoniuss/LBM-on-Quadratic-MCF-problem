function [mu, exitflag, theta, rho] = DSBM(B_z, B_alpha, l, mu)

muSize = size(mu,1);
muVar = sdpvar(muSize,1);
vVar = sdpvar(1,1);

Constraints = [ (l >= vVar):'rho', (vVar >= B_z'*muVar - B_alpha'):'theta'];
Objective = vVar + 0.5*norm(muVar - mu)^2;

%options = sdpsettings('verbose', 0, 'solver','quadprog');
options = sdpsettings('verbose', 0, 'solver','mosek');
%options = sdpsettings('verbose', 0);

sol = optimize(Constraints,Objective, options);
mu = value(muVar);
theta = dual(Constraint('theta'));
rho = dual(Constraint('rho'));
exitflag = sol.problem;

end
