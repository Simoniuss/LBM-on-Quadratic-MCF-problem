function [mu, exitflag, theta] = LBM_base(B_z, B_alpha, l, mu)

muSize = size(mu,1);

muVar = sdpvar(muSize,1);

%Constraints = B_z'*muVar - B_alpha' <= l;
%Constraints = [l * ones( size( B_z , 2 ) , 1 ) >=  B_z' * ( muVar) - B_alpha'];

%try to rewrite it
Constraints = [l  >= max(B_z' * ( muVar) - B_alpha')];


Objective = 0.5*norm(muVar - mu)^2;

%options = sdpsettings('verbose', 0, 'solver','quadprog');
options = sdpsettings('verbose', 0, 'solver','mosek');
%options = sdpsettings('verbose', 0);

sol = optimize(Constraints, Objective, options);
mu = value(muVar);
theta = dual(Constraints);
exitflag = sol.problem;

end
