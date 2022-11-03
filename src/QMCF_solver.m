function [x_sol] = QMCF_solver(Q, q, E, b, u, epsilon, l, lambda, best_l, m_lbm)
   [m,n] = size(E);
   [U,Sigma,V,U_m] = compactSVD(E);
   %get the full Q through Q = diag(Q)
   %tutti i vettori sono presi come riga
   if (b*U_m == 0)
       disp('Ex=b unsatisfiable')
       x_sol = zeros(n,1);
       return 
   end
   x_sat_const_found = false;
   bar_mu = zeros(m,1); 
   bar_x = getBoxedx(Q, q, E, b, u, bar_mu) % nx1
   x_best = zeros(n,1);
   dual_f = @Dual_f
   z = (E*bar_x)-b; %1xm
   alpha = (z*bar_mu) - dual_f(bar_x, Q, bar_mu, E, b); %scalar
   B_z = [z]; % matix of vectors, mx1 ( mxi, with i num of iterations)
   B_alpha = [alpha]; %vector of scalars
   X = []; %matrix of vectors, nx1
   x_lin_const = V * pinv(Sigma) * U'*b;

   num_iterations = 0;
   
   while(num_iterations < 1000 && x_sat_const_found == false)
       [bar_x, bar_mu, x_sat_const_found] = main_structure(Q, q, E, b, u, epsilon, l, lambda, best_l, m_lbm, bar_mu, X); 
       X = [X; bar_x];
       num_iterations = num_iterations +1;
   end

    if(x_sat_const_found == false)
        disp('x feasible not found, max num of iteration reached')
        return;
   else
        disp('x feasible found')
   end
    
   x_best = bar_x; %set the new x feasbile found, so the best so far

   LB = L (bar_x, Q, q, bar_mu, E, b);
   UB = f (bar_x, Q, q);

    
   %reset the value for the new loop
   x_sat_const_found == false;
   num_iterations = 0;

   while(abs(UB-LB) >= epsilon && num_iterations < 1000)
       [bar_x, bar_mu, x_sat_const_found] = main_structure(Q, q, E, b, u, epsilon, l, lambda, best_l, m_lbm, bar_mu, X);
       
       if(x_sat_const_found==true)
        [UB, LB, l, x_best] = update_values(bar_x, x_best);
        x_sat_const_found=false;
       end
       
   end
end

%%other tool functions


function [UB, LB, l, x_best] = update_values(bar_x, x_best)
    if(f (bar_x, Q, q) < f (x_best, Q, q)) %update 
        x_best = bar_x;
        l = f (x_best, Q, q);
        UB = f (x_best, Q, q);
    end

    if(L (bar_x, Q, q, bar_mu, E, b) > LB)
        LB = L (bar_x, Q, q, bar_mu, E, b);
    end
end

%implement instruction from 21 to 34, WHEN CALLED IN FIRST LOOP

function [bar_x, bar_mu, x_sat_const_found] = main_structure(Q, q, E, b, u, epsilon, l, lambda, best_l, m_lbm, bar_mu, X)
        x_sat_const_found = false;

        [bar_mu, theta] = LBM(bar_mu); %theta is a vector of scalars
        sum_zb_thetab = B_z * theta'; %compute z_b*theta_b -> out is a vec
        if(sum(theta)==1 && sum(sum_zb_thetab)==0)
          bar_x = X * theta';
          x_sat_const_found = true;    
        else %try with lagrangian heuristic
            if(sum(theta)~=1)
                theta = theta/sum(theta);
            end
            bar_x = X * theta';
            bar_x = bar_x - V*V'*bar_x - x_lin_const;

            x_sat_const_found = true; %%pre set to true, then check
            for i = 1:n
                if(bar_x(i)<0 || bar_x(i) > u)
                    x_sat_const_found = false;
                end
            end
            
            %%se x_sat_const_found false
            if(x_sat_const_found == false)
                bar_x = getBoxedx(Q, q, E, b, u, bar_mu);
            end
        end 

end



function [bar_mu, theta] = LBM(n)
    bar_mu = zeros(m,1)
    theta = ones(n,1)
end

function [x] = getBoxedx(Q, q, E, b, u, mu)
    % Return a feasible x which solves the the problem:
    %       min x'*Q*x + q*x + mu*( E*x - b )
    % Then project that x into the set [0,u]
    
    n = size(Q);
    x = zeros(n);
    ET = E';
    for i = 1:n
        if Q(i) == 0
            if mu'*b < q(i)*u(i) + mu'*( E(:,i)'*u(i) - b)
                x(i) = 0;
            else
                x(i) = u(i);
            end
        else
            x(i) = (-q(i) - ET(i,:)*mu)/(2*Q(i));
            x(i) = max(0, min(x(i), u(i)));
        end
    end
end

%parametric fun in mu
function [dual_f] = Dual_f(x, Q, mu, E, b)
    dual_f  = -x'*Q*x - q*x - mu'*(E*x - b) 
end                                          % la parentesi va trasposta per moltiplicarla con mu

%parametric fun in x and mu
function [L_x_mu] = L (x, Q, q, mu, E, b)
    L_x_mu = x'*Q*x + q*x + mu'*(E*x-b)
end


%parametric fun in x
function [f_x] = f (x, Q, q)
    f_x = x'*Q*x + q*x;
    return
end

