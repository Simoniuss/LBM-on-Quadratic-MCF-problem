function [x, exitFlag] = QMCF_solver_last_triple(Q, q, E, b, u, epsilon, l, lambda, best_l, m1, max_iter)

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
z = b - E*x;
f_val = -L(x,mu);
B_mu = mu;
B_z = z; 
B_f = f_val;

iter = 0;
SS = 0;
NS = 0;


while(iter < max_iter)
    
    % Compute new mu
    unbounded = true;
    while(unbounded)
        [d, exitflag, theta, f_B] = LBM_base_triple(B_z, B_f, B_mu, l, mu);
       % [newmu, exitflag, theta, f_B] = LBM_base_triple(B_z, B_f, B_mu, l, mu);
        %[newmu, exitflag, theta, rho] = DSBM(B_z, B_alpha, l, mu);
        
        if exitflag == 2
            disp('Unbounded')
            l = lambda*-L(x, mu) + (1 - lambda)*v;
        elseif exitflag == 0
            unbounded = false;
        else
            disp(exitflag)
            disp('Error')
            return
        end
    end
    
    newmu = mu+d;
    newx = getBoxedx(Q, q, E, b, u, newmu);

    disp("dist")
    disp(norm(newmu-mu))

    z = b - E*newx;
    f_val = (-L(newx, newmu));
    B_z = [ B_z z ];
    B_f = [ B_f f_val ];
    B_mu = [B_mu newmu];
    
    % mu update
    v = f_B;
    disp("V")
    disp(v)
    disp("Lnew")
    disp(f_val)
    if -L(newx, newmu) <= -L(x, mu) + m1*(v - (-L(x, mu)))
        step = 'SS';
        mu = newmu;
        x = newx;        
        % LS
        l = lambda*-L(x, mu) + (1 - lambda)*v;
        disp("NEWL")
        disp(l)
    else
        step = 'NS';
    end
       

    iter = iter + 1;
    if (step == "NS")
        NS = NS +1;
    else
        SS = SS +1;
    end
    
    fprintf('iter \t f-L \t\t norm(E*x-b) \t f \t\t\t\t\t L \t\t\t level \t\t\t step \t #bundle \t SS \t NS\n');
    fprintf('%d \t %1.4e \t %1.4e \t %1.4e \t %1.4e \t %1.4e \t %s \t %d \t\t\t %d \t\t %d\n', ...
        iter, f(x)-L(x,mu), norm(E*x-b), f(x), -L(x,mu), l, step, size(B_z,2), SS, NS);
end

if (iter == max_iter)
    disp('Maximum iterations exceeded');
    exitFlag = 2;
    return
end
end
