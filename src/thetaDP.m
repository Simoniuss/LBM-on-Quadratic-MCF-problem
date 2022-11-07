function [theta, exitflag] = thetaDP(B_z, B_alpha, l)
% Solve the DP for the LBM finding the weights theta used to build a
% a feasible solution. This function solves the following problem:
%       theta = arg min { 0.5*norm(sum(z^b*theta^b))^2 + sum((l+alpha^b)*theta), b \in B}
% B is the bundle which contains subgradients z and alpha, and l is the
% level parameter.
%

% Number of pairs in the bundle
bundleSize = size(B_z,2);

% theta
thetaVar = sdpvar(bundleSize,1);

Constraints = thetaVar >= 0;

Objective = 0.5*norm(B_z*thetaVar)^2 + (l+B_alpha)*thetaVar;

%options = sdpsettings('verbose', 0, 'solver','quadprog');
options = sdpsettings('verbose', 0, 'solver','mosek');
%options = sdpsettings('verbose', 0);

sol = optimize(Constraints,Objective, options);
theta = value(thetaVar);
exitflag = sol.problem;

end
