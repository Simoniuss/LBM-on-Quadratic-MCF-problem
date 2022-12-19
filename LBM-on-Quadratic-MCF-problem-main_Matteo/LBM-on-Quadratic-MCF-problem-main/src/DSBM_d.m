function [d, exitflag, theta] = DSBM_d(B_z, B_alpha, l, mu)

% Number of pairs in the bundle
bundleSize = size(B_z,2);
dSize = size(mu,1);

dVar = sdpvar(dSize,1);
vVar = sdpvar(1,1);


Constraints = [ l >= vVar, (vVar * ones(bundleSize,1) >= B_z'*(mu+dVar) - B_alpha'): 'bundle'];

Objective = vVar + 0.1*0.5*norm(dVar)^2;

%options = sdpsettings('verbose', 0, 'solver','quadprog');
options = sdpsettings('verbose', 0, 'solver','mosek');
%options = sdpsettings('verbose', 0);

%if bundleSize == 1
%    savecplexlp(Constraints, Objective, 'DSBM1');
%elseif bundleSize == 2
%    savecplexlp(Constraints, Objective, 'thetaDPData2');
%end

sol = optimize(Constraints,Objective, options);
d = value(dVar);
theta = dual(Constraints('bundle'));
exitflag = sol.problem;

end
