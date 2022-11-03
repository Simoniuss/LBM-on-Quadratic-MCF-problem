function [x_best, exitFlag] = QMCF_solver(Q, q, E, b, u, epsilon, l, lambda, best_l, m_lbm, max_iter)
   exitFlag=1;
   [m,n] = size(E);
   [U,Sigma,V,U_m] = compactSVD(E);
   %get the full Q through Q = diag(Q)
   %tutti i vettori sono presi come riga
   if (U_m'*b > 1.0000e-10)
       disp('Ex=b unsatisfiable');
       exitFlag = -1;
       return 
   end
   x_sat_const_found = false;
   bar_mu = zeros(m,1); 
   bar_x = getBoxedx(Q, q, E, b, u, bar_mu); % nx1
   x_best = [];
   %dual_f = @Dual_f;
   dualf = @(mu) -bar_x' * diag(Q) * bar_x - q' * bar_x - mu'*(E*bar_x - b);
   z = (E*bar_x)-b; %1xm
   %alpha = dot(z, bar_mu) - dual_f(bar_x, Q, q, bar_mu, E, b); %scalar
   alpha = dot(z, bar_mu) - dualf(bar_mu);
   B_z = [z]; % matix of vectors, mx1 ( mxi, with i num of iterations)
   B_alpha = [alpha]; %vector of scalars
   X = [bar_x]; %matrix of vectors, nx1
   x_lin_const = V * pinv(Sigma) * U'*b;
   num_iterations = 0;    
   x_best = zeros(n,1); %set the new x feasbile found, so the best so far
   
   LB = 0;
   UB = 2*epsilon;
   %LB = L (bar_x, Q, q, bar_mu, E, b);
   %UB = f (bar_x, Q, q);

   
   %reset the value for the new loop
   x_sat_const_found = false;
   num_iterations = 0;

   while(abs(UB-LB) >= epsilon && num_iterations < max_iter)
        
        [bar_mu, theta, l] = LBM(dualf, bar_mu, l, B_z, B_alpha, lambda, best_l, m_lbm);
        %[bar_mu, theta] = LBM(bar_mu); %theta is a vector of scalars
        if(sum(theta)==1 && ~any(B_z * theta)) % Case 1
          bar_x = X * theta;
          x_sat_const_found = true;    
        else %try with lagrangian heuristic (Case 2)
            if(sum(theta)~=1 && sum(theta)~=0)
                theta = theta/sum(theta);
            end
            bar_x = X * theta;
            bar_x = bar_x - V*V'*bar_x - x_lin_const;

            if (all(bar_x >= zeros(n,1)) && all(bar_x <= u))
                x_sat_const_found = true;
            else
                bar_x = getBoxedx(Q, q, E, b, u, bar_mu);
            end
        end 
       
       if(x_sat_const_found)
           if isempty(x_best)
               x_best = bar_x;
           else
               if (f (bar_x, Q, q) < f (x_best, Q, q)) %update
                    x_best = bar_x;
                    l = f (x_best, Q, q);
                    UB = l;
               end
           end
           if(L (bar_x, Q, q, bar_mu, E, b) > LB)
                LB = L (bar_x, Q, q, bar_mu, E, b);
           end
           x_sat_const_found=false;
       end
       
        X = [X bar_x];
        % Compute new pair (z, alpha) and append to the bundle
        newz = B_z*theta;
        newalpha = B_alpha*theta;
        B_z = [ B_z newz ];
        B_alpha = [ B_alpha newalpha ];
        %disp([num_iterations, theta]);
        num_iterations = num_iterations +1;
       
   end

   if (num_iterations == max_iter)
       disp('Maximum iterations exceeded');
       exitFlag = 2;
       return 
   end
end
   
%%other tool functions


%parametric fun in x and mu
function [L_x_mu] = L (x, Q, q, mu, E, b)
    L_x_mu = x'*diag(Q)*x + q'*x + mu'*(E*x-b);
end


%parametric fun in x
function [f_x] = f (x, Q, q)
    f_x = x'*diag(Q)*x + q'*x;
    return
end

