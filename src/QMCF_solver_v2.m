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

% Initialization of the bundle
dualf = @(barx, mu) -L(barx, mu); %inv fun via - sign
z = b - (E*bar_x); %mx1
dualf_mu = dualf(bar_x,bar_mu);
alpha = dot(z, bar_mu) - dualf_mu;
B_z = z; % matrix of vectors, mx1 ( mxi, with i num of iterations)
B_alpha = alpha; %vector of scalars



%Initialization of variable for optimal x
num_iterations = 0;


fprintf('iter \t ||bar_mu|| \t ||bar_x|| \t l \t ||New z|| \t New alpha\n');

while(num_iterations < max_iter)
    
    if dualf_mu < l
        l = dualf_mu-l;
    end
    
    unbounded = true;
    while(unbounded)
        [bar_mu, exitflag] = LBM_v2(B_z, B_alpha, l, bar_mu);
        if exitflag == 2
            l = lambda*dualf(bar_x,bar_mu) - (1 - lambda)*l;
        elseif exitflag == 0
            unbounded = false;
        else
            disp('Error')
            return
        end
    end

    bar_x = getBoxedx(Q, q, E, b, u, bar_mu);
    
    %disp(f(bar_x))
    %disp(dualf(bar_x,bar_mu))

    newz = -(E*bar_x)+b;
    dualf_mu = dualf(bar_x,bar_mu);
    newalpha = newz'*bar_mu - dualf_mu;
    
    B_z = [ B_z newz ];
    B_alpha = [ B_alpha newalpha ];
    num_iterations = num_iterations +1;
    
    disp(bar_x)
    
    
    fprintf('%d \t %1.4e \t %1.4e \t %1.4e \t %1.4e \t %1.4e\n', ...
        num_iterations, norm(bar_mu), norm(bar_x), l, norm(newz), newalpha);
    
end

if (num_iterations == max_iter)
    disp('Maximum iterations exceeded');
    exitFlag = 2;
    x_best=bar_x;
    return
end
end

