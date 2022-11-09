function [x_best, exitFlag] = QMCF_solver(Q, q, E, b, u, epsilon, l, lambda, best_l, m_lbm, max_iter)

% exitflag: 1: optimal found, -1: unsatisfiable, 2: max iter reached
exitFlag=1;
[m,n] = size(E);
[U,Sigma,V,U_m] = compactSVD(E);

% Check satisfiability condition for E*x=b
if (U_m'*b > 1.0e-10)
    disp('Ex=b unsatisfiable');
    exitFlag = -1;
    return
end

% Initilization of parameters
f = @(x) x'*diag(Q)*x + q'*x;
bar_mu = zeros(m,1);
bar_x = getBoxedx(Q, q, E, b, u, bar_mu); % nx1
x_best = [];
L = @(x, mu) f(x) + mu'*(E*x-b);
x_lin_const = V * pinv(Sigma) * U'*b;

% Initialization of the bundle
dualf = @(barx, mu) -barx' * diag(Q) * barx - q' * barx - mu'*(E*barx - b);
z = (E*bar_x)-b; %mx1
alpha = dot(z, bar_mu) - dualf(bar_x, bar_mu);
B_z = z; % matrix of vectors, mx1 ( mxi, with i num of iterations)
B_alpha = alpha; %vector of scalars
X = bar_x; %matrix of vectors, nx1

% Initialization of lower bound and upper bound
LB = 0;
UB = Inf;
%LB = L(bar_x, bar_mu);
%UB = f(bar_x);


%Initialization of variable for optimal x
x_sat_const_found = false;
num_iterations = 0;
theta_case = 0;
normd = 0;

fprintf('iter \t gap \t ||theta|| \t ||bar_mu|| \t ||bar_x|| \t l \t Case \t ||New z|| \t New alpha \t ||d|| \t step\n');

while(abs(UB-LB) >= epsilon && num_iterations < max_iter)
    
    [bar_mu, theta, l, normd, steptype] = LBM(dualf, bar_x, bar_mu, l, B_z, B_alpha, lambda, best_l, m_lbm);
    
    % Case 1: theta is a convex combinator, optimal x
    if(sum(theta)==1 && all((B_z * theta)<eps)) % Case 1
        theta_case = 1;
        bar_x = X * theta;
        x_sat_const_found = true;
    
    % Case 2: try to force x in [0,u] using a Lagrangian heuristic
    else
        theta_scaled = theta;
        if(sum(theta)~=1 && sum(theta)~=0)
            theta_scaled = theta/sum(theta);
        end
        
        bar_x = X * theta_scaled;
        bar_x = bar_x - V*V'*bar_x - x_lin_const;
        
        disp(x_lin_const')
        
        if (all(bar_x >= 0) && all(bar_x <= u))
            theta_case = 2;
            theta = theta_scaled;
            x_sat_const_found = true;
        
        % Case 3: optimal x not found, find a feasible one
        else
            theta_case = 3;
            bar_x = getBoxedx(Q, q, E, b, u, bar_mu);
        end
    end
    
    
    % Update LB and UB if an optimal x has been found
    if(x_sat_const_found)
        if isempty(x_best)
            x_best = bar_x;
        else
            
            % Update UB and the level parameter
            if (f(bar_x) < f(x_best))
                x_best = bar_x;
                l = f (x_best);
                UB = l;
            end
        end
        
        % Update LB
        if(L(bar_x, bar_mu) > LB)
            LB = L(bar_x, bar_mu);
        end
        x_sat_const_found=false;
    end
    
    X = [X bar_x];
    % Compute new pair (z, alpha) and append to the bundle
    newz = B_z*theta;
    newalpha = B_alpha*theta;
    B_z = [ B_z newz ];
    B_alpha = [ B_alpha newalpha ];
    num_iterations = num_iterations +1;
    
    
    fprintf('%d \t %1.1e \t %1.1f \t %1.1e \t %1.1e \t %1.1e \t %d \t %1.1e \t %1.1e \t %1.1e \t %s\n', ...
        num_iterations, abs(UB-LB), norm(theta), norm(bar_mu), norm(bar_x), l, theta_case, norm(newz), newalpha, normd, steptype);
    
end

if (num_iterations == max_iter)
    disp('Maximum iterations exceeded');
    exitFlag = 2;
    return
end
end

