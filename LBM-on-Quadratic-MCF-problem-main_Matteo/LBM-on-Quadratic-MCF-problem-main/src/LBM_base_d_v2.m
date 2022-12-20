function [d, exitflag, theta] = LBM_base_d_v2(B_z, B_alpha, l, mu)

dSize = size(mu,1);

dVar = sdpvar(dSize,1);

Constraints = [];

%const of FRNAGIONI PBM version
%Constraints = [l * ones( size( B_z , 2 ) , 1 ) >=  B_alpha' + B_z' * (mu+dVar)];

%OK working, done with SIMO
Constraints = B_z'*(mu+dVar) - B_alpha' <= l;

%rewrite of previous const done with SIMO
%Constraints = [l * ones( size( B_z , 2 ) , 1 ) >=  B_z' * ( dVar) - B_alpha'];

%try to rewrite it again
%Constraints = [l  >= max(B_z' * ( mu+dVar) - B_alpha')];

Objective = norm(dVar)^2 / 2;

%options = sdpsettings('verbose', 0, 'solver','quadprog');
options = sdpsettings('verbose', 0, 'solver','mosek');
%options = sdpsettings('verbose', 0);

%SI PUO ESTRARRE V DA QUA ?
sol = optimize(Constraints, Objective, options);
d = value(dVar);
theta = dual(Constraints);
exitflag = sol.problem;

end
