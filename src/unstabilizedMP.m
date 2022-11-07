function [mu, v, exitflag] = unstabilizedMP(B_z, B_alpha)
% Solve the unstabilezed MP for the LBM. It is a LP problem in the 
% following form:
%       (mu, v) = arg min {v : v >= dot(z^b, mu) - alpha^b, b \in B}
% B is the bundle which contains subgradients z and alpha.
%

% Number of pairs in the bundle
bundleSize = size(B_z,2); %num of columns
muSize = size(B_z,1);

muVar = sdpvar(muSize,1);
vVar = sdpvar(1,1);

Constraints = [];
for i = 1:bundleSize
    Constraints = [Constraints, vVar >= B_z(:, i)'*muVar - B_alpha(i)];
end

Objective = vVar;

%options = sdpsettings('verbose', 0, 'solver','linprog');
options = sdpsettings('verbose', 0, 'solver','mosek');
%options = sdpsettings('verbose', 0);

sol = optimize(Constraints,Objective, options);

mu = value(muVar); %make the variables printable
v = value(vVar);
exitflag = sol.problem;

end

