function [d, exitflag] = TR_d(B_z, B_alpha, delta, mu)

dSize = size(mu,1);

% mu
dVar = sdpvar(dSize,1);
vVar = sdpvar(1,1);

Constraints = [vVar >= B_z'*(mu+dVar) - B_alpha', norm(dVar) <= delta];

Objective = vVar;

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
exitflag = sol.problem;

end
