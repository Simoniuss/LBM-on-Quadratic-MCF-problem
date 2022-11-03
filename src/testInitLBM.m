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
    %def variables
mu = optimvar('mu', size(b));
v = optimvar('v');
    %solve problem
prob = optimproblem;
prob.Objective = v;
    %define the # of constraints (= dim del bundle)
bundleconstr = optimconstr(size(B_z,2));
    %define the actual constraints
for i = 1:size(B_z,2)
    bundleconstr(i) = v >= B_z(:, i)'*mu - B_alpha(i);
end
    %put constraints defined in the variable for the optimier
prob.Constraints.cons = bundleconstr;

% if unbounded exitflag = -3
[sol, fval, exitflag] = solve(prob); 

%sol = (mu,v)
%fval = func val, equal to v


l = 0; %init since can be unbounded, so it could have no value

niter = size(B_z,2);

% theta problem
theta = optimvar('theta', niter);
prob = optimproblem;
prob.Objective = sum((B_z*theta).^2) + ((l+B_alpha)*theta);%here opt find theta
prob.Constraints.cons = theta >= 0;
% if unbounded exitflag = -3
[sol, fval, exitflag] = solve(prob);

theta = sol.theta;

d = - B_z*theta;
v = max(B_z'*(barmu + d) - B_alpha);

%condition for NS or SS
if dualf(barx, barmu+d) - dualf(barx, barmu) <= 0.001*(v- dualf(barx, barmu))
    barmu = barmu+d;
end

newz = B_z*theta;
newalpha = B_alpha*theta;

B_z = [B_z newz];
B_alpha = [B_alpha newalpha];





