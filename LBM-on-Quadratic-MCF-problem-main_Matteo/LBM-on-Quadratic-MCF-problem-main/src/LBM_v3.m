function [d, exitflag, theta] = LBM_v3(B_z, B_alpha, l, mu)

% Number of pairs in the bundle
dSize = size(mu,1);

% mu
dVar = sdpvar(dSize,1);

Constraints = [B_z'*(mu+dVar) - B_alpha' <= l];

Objective = 0.5*norm(dVar)^2;

%options = sdpsettings('verbose', 0, 'solver','quadprog');
options = sdpsettings('verbose', 0, 'solver','mosek');
%options = sdpsettings('verbose', 0);

%if bundleSize == 1
%savecplexlp(Constraints, Objective, 'LBM_v2_1');
%elseif bundleSize == 2
%    savecplexlp(Constraints, Objective, 'thetaDPData2');
%end

sol = optimize(Constraints,Objective, options);
d = value(dVar);
theta = dual(Constraints);
exitflag = sol.problem;

end
