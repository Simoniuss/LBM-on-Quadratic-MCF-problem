testScript;
% Initialization
dualf = @(x, mu) -x' * diag(Q) * x - q' * x - mu'*(E*x - b);
barmu = zeros(size(b));
barx = getBoxedx(Q, q, E, b, u, barmu);
z = E*barx - b;
alpha = dot(z, barmu) - dualf(barx, barmu);
B_z = z ; 
B_alpha = alpha ;

% bestl problem
mu = optimvar('mu', size(b));
v = optimvar('v');
prob = optimproblem;
prob.Objective = v;
bundleconstr = optimconstr(size(B_z,2));
for i = 1:size(B_z,2)
    bundleconstr(i) = v >= B_z(:, i)'*mu - B_alpha(i);
end
prob.Constraints.cons = bundleconstr;
% if unbounded exitflag = -3
[sol, fval, exitflag] = solve(prob);

l = 0;

niter = size(B_z,2);

% theta
theta = optimvar('theta', niter);
prob = optimproblem;
prob.Objective = sum(sum(B_z*theta, 2).^2) + sum((l+B_alpha)*theta);
prob.Constraints.cons = theta >= 0;
% if unbounded exitflag = -3
[sol, fval, exitflag] = solve(prob);

theta = sol.theta;

d = - sum(B_z*theta, 2);
v = max(B_z'*(barmu + d) - B_alpha);

if dualf(barx, barmu+d) - dualf(barx, barmu) <= 0.001*(v- dualf(barx, barmu))
    barmu = barmu+d;
end

newz = sum(B_z*theta, 2);
newalpha = sum(B_alpha*theta);

B_z = [B_z newz];
B_alpha = [B_alpha newalpha];





