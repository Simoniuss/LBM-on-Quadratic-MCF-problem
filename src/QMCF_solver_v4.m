function [x_best, exitFlag] = QMCF_solver_v4(Q, q, E, b, u, epsilon, l, lambda, best_l, m_lbm, max_iter)

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

% Best parameters
best_x = [1; 2.838; 6.162; 0; 1.0; 2.838];
best_fx = 6.3282e+03;
%best_fx = l;
%best_x_founded = [1.0001; 2.8380; 6.1619; 0; 1.0002; 2.8379];
%best_fx_founded = 6.6455e+03;

% Initilization of parameters
f = @(x) x'*diag(Q)*x + q'*x;
bar_mu = zeros(m,1);
bar_x = getBoxedx(Q, q, E, b, u, bar_mu); % nx1
L = @(x, mu) f(x) + mu'*(E*x-b);
x_lin_const = V * pinv(Sigma) * U'*b;
x_best = [];

% Initialization of the bundle
dualf = @(barx, mu) -L(barx, mu); %inv fun via - sign
z = b - (E*bar_x); %mx1
dualf_mu = dualf(bar_x,bar_mu);
alpha = dot(z, bar_mu) - dualf_mu;
B_z = z; % matrix of vectors, mx1 ( mxi, with i num of iterations)
B_alpha = alpha; %vector of scalars
X = bar_x; %matrix of vectors, nx1

% Initialization of lower bound and upper bound
LB = -Inf;
UB = Inf;


%Initialization of variable for optimal x
x_sat_const_found = false;
num_iterations = 0;
theta_case = 0;

f_L_gap = Inf;

%close
%tiledlayout(3,1)

fprintf('iter \t UB-LB \t norm(E*x-b) \t f \t\t L \t\t Case \t #bundle \t Gap\n');

while(abs(UB-LB) >= epsilon && num_iterations < max_iter)
    
    % LS
    %l = lambda*dualf(bar_x, bar_mu) + (1 - lambda)*-best_fx;
    l = lambda*dualf(bar_x, bar_mu) + (1 - lambda)*l;
    %l = -best_fx;
    
    %tic;
    % Compute new mu
    unbounded = true;
    while(unbounded)
        %[d, exitflag, theta] = LBM_v3(B_z, B_alpha, l, bar_mu);
        %[d, exitflag] = TR_d(B_z, B_alpha, delta, bar_mu);
        [d, exitflag, theta] = DSBM_d(B_z, B_alpha, l, bar_mu);
        
        if exitflag == 2
            l = lambda*dualf(bar_x, bar_mu) + (1 - lambda)*-best_fx;
        elseif exitflag == 0
            unbounded = false;
        else
            disp(exitflag)
            disp('Error')
            return
        end
    end
    %time = toc;
    
    %fprintf('Yalmip: %1.8e\n', time);
    
    %tic;
    % Add new point to the bundle
    newz = b - (E*bar_x);
    dualf_mu = dualf(bar_x, bar_mu+d);
    newalpha = newz'*(bar_mu+d) - dualf_mu;
    
    % Compute SS or NS on mu
    v = max(B_z'*(bar_mu+d) - B_alpha');
    if dualf(bar_x, bar_mu+d) < dualf(bar_x, bar_mu) + m_lbm*(v - dualf(bar_x, bar_mu))
        steptype = 'SS';
        bar_mu = bar_mu+d;
    else
        steptype = 'NS';
    end
    
    % Case 1: theta is a convex combinator, optimal x
%     disp('Case 1 check')
%     disp(sum(theta))
%     disp(norm(B_z * theta))

    if(abs(sum(theta)-1)<=epsilon && all(abs(B_z * theta)<=epsilon)) % Case 1
        theta_case = 1;
        new_bar_x = X * theta;
        x_sat_const_found = true;
    
    % Case 2: try to force x in [0,u] using a Lagrangian heuristic
    else
        theta_scaled = theta;
        if(sum(theta)~=1 && sum(theta)~=0)
            theta_scaled = theta/sum(theta);
        end
        
        new_bar_x = X * theta_scaled;
        new_bar_x = new_bar_x - V*V'*new_bar_x + x_lin_const;
        
        
        if (all(new_bar_x >= -epsilon) && all(new_bar_x <= u+epsilon))
            theta_case = 2;
            theta = theta_scaled;
            x_sat_const_found = true;
        
        % Case 3: optimal x not found, find a feasible one
        else
            theta_case = 3;
            new_bar_x = getBoxedx(Q, q, E, b, u, bar_mu);
        end
    end
    
    
    %if f(new_bar_x) > f(bar_x)
    if abs(bar_mu'*(E*new_bar_x-b)) < f_L_gap
        bar_x = new_bar_x;
        % Bundle reset
        newz = b - (E*bar_x);
        dualf_mu = dualf(bar_x, bar_mu);
        newalpha = newz'*bar_mu - dualf_mu;
        B_z = newz;
        B_alpha = newalpha;
        X = bar_x;
        f_L_gap = abs(bar_mu'*(E*bar_x-b));
    else
        % Bundle update 
        B_z = [ B_z newz ];
        B_alpha = [ B_alpha newalpha ];
        X = [X new_bar_x];
    end
     
    % Update LB and UB if an optimal x has been found
    if(x_sat_const_found)
        if isempty(x_best)
            x_best = bar_x;
            LB = L(x_best, bar_mu);
            UB = f(x_best);
        else
            
            % Update UB and the level parameter
            if (f(bar_x) < f(x_best))
                x_best = bar_x;
                UB = f(x_best);
            end
            
            % Update LB
            if(L(bar_x, bar_mu) > LB)
                LB = L(bar_x, bar_mu);
            end
        end
        
        x_sat_const_found=false;       
    end
    
    
    % Compute new pair (z, alpha) and append to the bundle
%     if(abs(sum(theta)-1)<=epsilon)
%         newz = B_z*theta;
%         newalpha = B_alpha*theta;
%     else
%         newz = -(E*bar_x)+b;
%         newalpha = newz'*bar_mu - dualf(bar_x, bar_mu);
%     end
    
    
    
    num_iterations = num_iterations + 1;
    %time = toc;
    
    %fprintf('Algorithm: %1.8e\n', time);
    
    
    fprintf('%d \t %1.4e \t %1.4e \t %1.4e \t %1.4e \t %d \t %d \t %1.4e\n', ...
        num_iterations, abs(UB-LB), norm(E*bar_x-b), f(bar_x), L(bar_x,bar_mu), theta_case, size(B_z,2), bar_mu'*(E*bar_x-b));
    
    % print informations
    %fprintf('iter \t UB-LB \t\t norm(E*x-b) \t ||bar_mu|| \t ||bar_x|| \t l \t\t dualf \t\t step \t ||New z|| \t New alpha\n');
    %fprintf('%d \t %1.4e \t %1.4e \t %1.4e \t %1.4e \t %1.4e \t %1.4e \t %s \t %1.4e \t %1.4e\n\n', ...
    %    num_iterations, abs(UB-LB), norm(E*bar_x-b), norm(bar_mu), norm(bar_x), l, dualf_mu, steptype, norm(newz), newalpha);
    
%     tic;
%     % Plots
%     nexttile(1)
%     title('|| E*x-b ||')
%     plot(num_iterations, norm(E*bar_x-b), 'o')
%     hold on
%     nexttile(2)
%     title('Relative gap')
%     plot(num_iterations, abs(f(bar_x) - best_fx)/best_fx, 'o')
%     hold on
% %     nexttile(2)
% %     title('|| x - x* ||')
% %     plot(num_iterations, norm(bar_x - best_x), 'o')
% %     hold on
%     nexttile(3)
%     title('|UB - LB|')
%     plot(num_iterations, abs(UB-LB), 'o')
%     hold on
%     time = toc;
%     
%     fprintf('Plot: %1.8e\n', time);
    
end
%hold off

if (num_iterations == max_iter)
    disp('Maximum iterations exceeded');
    exitFlag = 2;
    return
end
end

