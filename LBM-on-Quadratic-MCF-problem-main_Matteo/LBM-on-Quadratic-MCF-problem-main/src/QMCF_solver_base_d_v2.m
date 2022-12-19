function [x, exitFlag] = QMCF_solver_base_d_v2(Q, q, E, b, u, epsilon, l, lambda, best_l, m1, max_iter)

% exitflag: 1: optimal found, -1: unsatisfiable, 2: max iter reached
exitFlag=1;
[m,n] = size(E);
[U,S,V,um] = compactSVD(E);

% Check satisfiability condition for E*x=b
if (um'*b > 1.0e-10)
    disp('Ex=b unsatisfiable');
    exitFlag = -1;
    return
end

% Initilization of parameters
f = @(x) x'*diag(Q)*x + q'*x;
mu = zeros(m,1);
x = getBoxedx(Q, q, E, b, u, mu);
L = @(x, mu) f(x) + mu'*(E*x - b);

% Initialization of the bundle
z = b - E*x; % -L gradient
alpha = z'*mu - (-L(x,mu));
B_z = z; 
B_alpha = alpha;

iter = 0;
best_fx = -6.3282e+03;
l = lambda*-L(x, mu) + (1 - lambda)*l;

fprintf('iter \t f-L \t\t norm(E*x-b) \t f \t\t L \t\t level \t\t step \t #bundle\n');

while(iter < max_iter)
    
    % Compute new mu
    unbounded = true;
    while(unbounded)
        [d, exitflag, theta] = LBM_base_d_v2(B_z, B_alpha, l, mu, (-L(x,mu)));
        %[d, exitflag, v, theta, rho] = DSBM_base_d_v2(B_z, B_alpha, l, mu);
        
        if exitflag == 2
            disp('Unbounded')
            l = lambda*-L(x, mu) + (1 - lambda)*l;
        elseif exitflag == 0
            unbounded = false;
        else
            disp(exitflag)
            disp('Error')
            return
        end
    end
    
    % mu update
    v = max(B_z'*(mu+d) - B_alpha');
    if -L(x, mu+d) < -L(x, mu) + m1*(v - (-L(x, mu)))
        step = 'SS';
        mu = mu+d;
        
        % if SS compute the new x and check if it is better than the
        % previous one wrt Ex=b. If it's better reset the bundle
        xnew = getBoxedx(Q, q, E, b, u, mu);
        if norm(E*xnew - b) < norm(E*x - b)
            x = xnew;
            B_z = [];
            B_alpha = [];
        end
        
        z = b - E*x;
        alpha = z'*mu - (-L(x,mu));
        B_z = [ B_z z ];
        B_alpha = [ B_alpha alpha ];
        %B_alpha = (B_z'*(-d) + B_alpha(end) + (-L(x,mu)) - (-L(xprec,mu-d)))';
        
        % LS
        l = lambda*-L(x, mu) + (1 - lambda)*l;
        
    else
        % mu nell'iterazione successiva non cambia
        % però nel bundle ci metto mu+d
        step = 'NS';
        z = b - E*x;
        alpha = z'*(mu+d) - (-L(x,mu+d));
        B_z = [ B_z z ];
        B_alpha = [ B_alpha alpha ];
    end
    
       
    iter = iter + 1;
    
    
    fprintf('%d \t %1.4e \t %1.4e \t %1.4e \t %1.4e \t %1.4e \t %s \t %d\n', ...
        iter, f(x)-L(x,mu), norm(E*x-b), f(x), L(x,mu), l, step, size(B_z,2));
end

if (iter == max_iter)
    disp('Maximum iterations exceeded');
    exitFlag = 2;
    return
end
end

