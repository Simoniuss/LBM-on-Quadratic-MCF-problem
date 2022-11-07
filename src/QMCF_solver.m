function [x_best, exitFlag] = QMCF_solver(Q, q, E, b, u, epsilon, l, lambda, best_l, m, max_iter)
exitFlag=1;
[m,n] = size(E);
[U,Sigma,V,U_m] = compactSVD(E);
%get the full Q through Q = diag(Q)
%tutti i vettori sono presi come riga
if (U_m'*b > 1.0e-10)
    disp('Ex=b unsatisfiable');
    exitFlag = -1;
    return
end


f = @(x) x'*diag(Q)*x + q'*x;
x_sat_const_found = false;
bar_mu = zeros(m,1);
bar_x = getBoxedx(Q, q, E, b, u, bar_mu); % nx1
x_best = [];
L = @(x, mu) x'*diag(Q)*x + q'*x + mu'*(E*x-b);
dualf = @(mu) -bar_x' * diag(Q) * bar_x - q' * bar_x - mu'*(E*bar_x - b);
z = (E*bar_x)-b; %mx1
alpha = dot(z, bar_mu) - dualf(bar_mu);
B_z = [z]; % matrix of vectors, mx1 ( mxi, with i num of iterations)
B_alpha = [alpha]; %vector of scalars
X = [bar_x]; %matrix of vectors, nx1
x_lin_const = V * pinv(Sigma) * U'*b;
num_iterations = 0;

LB = 0;
UB = 2*epsilon;
%LB = L (bar_x, Q, q, bar_mu, E, b);
%UB = f (bar_x, Q, q);


%reset the value for the new loop
x_sat_const_found = false;
num_iterations = 0;

while(abs(UB-LB) >= epsilon && num_iterations < max_iter)
    
    [bar_mu, theta, l] = LBM(dualf, bar_mu, l, B_z, B_alpha, lambda, best_l, m);
    disp(['theta: ',num2str(nnz(theta))]);
    disp(['bar_mu: ',num2str(norm(bar_mu))]);
    %[bar_mu, theta] = LBM(bar_mu); %theta is a vector of scalars
    if(sum(theta)==1 && ~any(B_z * theta)) % Case 1
        disp('1');
        bar_x = X * theta;
        x_sat_const_found = true;
    else %try with lagrangian heuristic (Case 2)
        if(sum(theta)~=1 && sum(theta)~=0)
            theta = theta/sum(theta);
        end
        bar_x = X * theta;
        bar_x = bar_x - V*V'*bar_x - x_lin_const;
        
        if (all(bar_x >= 0) && all(bar_x <= u))
            disp('2');
            x_sat_const_found = true;
        else
            disp('3');
            bar_x = getBoxedx(Q, q, E, b, u, bar_mu);
        end
    end
    
    disp(['bar_x: ',num2str(nnz(bar_x))]);
    
    if(x_sat_const_found)
        if isempty(x_best)
            x_best = bar_x;
        else
            if (f(bar_x) < f(x_best)) %update
                x_best = bar_x;
                l = f (x_best);
                UB = l;
            end
        end
        if(L(bar_x, bar_mu) > LB)
            LB = L(bar_x, bar_mu);
        end
        x_sat_const_found=false;
    end
    
    dualf = @(mu) -bar_x' * diag(Q) * bar_x - q' * bar_x - mu'*(E*bar_x - b);
    X = [X bar_x];
    % Compute new pair (z, alpha) and append to the bundle
    newz = B_z*theta;
    newalpha = B_alpha*theta;
    B_z = [ B_z newz ];
    B_alpha = [ B_alpha newalpha ];
    num_iterations = num_iterations +1;
    
end

if (num_iterations == max_iter)
    disp('Maximum iterations exceeded');
    exitFlag = 2;
    return
end
end

