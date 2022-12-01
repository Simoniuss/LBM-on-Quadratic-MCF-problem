function [mu, exitflag] = DSBM(B_z, B_alpha, l, mu)
% Solve the DP for the LBM finding the weights theta used to build a
% a feasible solution. This function solves the following problem:
%       theta = arg min { 0.5*norm(sum(z^b*theta^b))^2 + sum((l+alpha^b)*theta), b \in B}
% B is the bundle which contains subgradients z and alpha, and l is the
% level parameter.
%

% Number of pairs in the bundle
bundleSize = size(B_z,2);
muSize = size(mu,1);


% mu
muVar = sdpvar(muSize,1);
vVar = sdpvar(1,1);


Constraints = [ l >= vVar, vVar * ones(bundleSize,1) >= B_z'*muVar - B_alpha'];

%Constraints = B_z'*(muVar) - B_alpha <= l;
Objective = vVar + 0.5*norm(muVar - mu)^2;

%options = sdpsettings('verbose', 0, 'solver','quadprog');
options = sdpsettings('verbose', 0, 'solver','mosek');
%options = sdpsettings('verbose', 0);

%if bundleSize == 1
%    savecplexlp(Constraints, Objective, 'DSBM1');
%elseif bundleSize == 2
%    savecplexlp(Constraints, Objective, 'thetaDPData2');
%end

sol = optimize(Constraints,Objective, options);
mu = value(muVar);
exitflag = sol.problem;

end
