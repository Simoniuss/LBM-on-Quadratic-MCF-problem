function [x_best, exitFlag] = QMCF_solver_v2(Q, q, E, b, u, epsilon, l, lambda, best_l, m_lbm, max_iter)

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

% Initilization of parameters
f = @(x) x'*diag(Q)*x + q'*x;
bar_mu = zeros(m,1);
bar_x = getBoxedx(Q, q, E, b, u, bar_mu); % nx1
L = @(x, mu) f(x) + mu'*(E*x-b);

% Initialization of the bundle
dualf = @(barx, mu) -L(barx, mu); %inv fun via - sign
z = b - (E*bar_x); %mx1
dualf_mu = dualf(bar_x,bar_mu);
alpha = dot(z, bar_mu) - dualf_mu;
B_z = z; % matrix of vectors, mx1 ( mxi, with i num of iterations)
B_alpha = alpha; %vector of scalars

% Lower bound and upper bound initialization
LB = -Inf;
UB = Inf;

%Initialization of flag variable
num_iterations = 0;

% Close figure and set new one
close
tiledlayout(3,1)

while(abs(UB-LB) >= epsilon && num_iterations < max_iter)
    
    % LS
    %l = lambda*dualf(bar_x, bar_mu) + (1 - lambda)*-best_fx;
    %l = lambda*dualf(bar_x, bar_mu) + (1 - lambda)*l;
    l = -best_fx;
    
    
    
        % TR example
    if num_iterations < 50
        d = 20;
    elseif num_iterations < 80
        d = 5;
    elseif num_iterations < 100
        d = 1;
    else
        d = 0.75;
    end

    % Compute new mu
    unbounded = true;
    while(unbounded)
        [new_bar_mu, exitflag] = LBM_v2(B_z, B_alpha, l, bar_mu);
        %[new_bar_mu, exitflag] = TR(B_z, B_alpha, d, bar_mu);
        %[new_bar_mu, exitflag] = DSBM(B_z, B_alpha, l, bar_mu);
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
    
    
    % Add new point to the bundle
    newz = b - (E*bar_x);
    dualf_mu = dualf(bar_x, new_bar_mu);
    newalpha = newz'*new_bar_mu - dualf_mu;
    
    % Compute SS or NS on mu
    v = max(B_z'*(new_bar_mu) - B_alpha');
    if dualf(bar_x, new_bar_mu) < dualf(bar_x, bar_mu) + m_lbm*(v - dualf(bar_x, bar_mu))
        steptype = 'SS';
        bar_mu = new_bar_mu;
    else
        steptype = 'NS';
    end
    
    % Compute SS or NS on x
    new_bar_x = getBoxedx(Q, q, E, b, u, bar_mu);
    
    if L(new_bar_x, bar_mu) < L(bar_x, bar_mu)
        bar_x = new_bar_x;
        % Bundle reset
        newz = b - (E*bar_x);
        dualf_mu = dualf(bar_x, bar_mu);
        newalpha = newz'*bar_mu - dualf_mu;
        B_z = [ newz ];
        B_alpha = [ newalpha ];
    else
        % Bundle update 
        B_z = [ B_z newz ];
        B_alpha = [ B_alpha newalpha ];
    end

     
    % LB and UB update
    LB = L(bar_x, bar_mu);
    UB = f(bar_x);

    
    num_iterations = num_iterations + 1;
    
    % print informations
    fprintf('iter \t UB-LB \t\t norm(E*x-b) \t ||bar_mu|| \t ||bar_x|| \t l \t\t dualf \t\t step \t ||New z|| \t New alpha\n');
    fprintf('%d \t %1.4e \t %1.4e \t %1.4e \t %1.4e \t %1.4e \t %1.4e \t %s \t %1.4e \t %1.4e\n\n', ...
        num_iterations, abs(UB-LB), norm(E*bar_x-b), norm(bar_mu), norm(bar_x), l, dualf_mu, steptype, norm(newz), newalpha);
    
    % Plots
    nexttile(1)
    title('|| E*x-b ||')
    plot(num_iterations, norm(E*bar_x-b), 'o')
    hold on
    nexttile(2)
    title('Relative gap')
    plot(num_iterations, abs(f(bar_x) - best_fx)/best_fx, 'o')
    hold on
%     nexttile(2)
%     title('|| x - x* ||')
%     plot(num_iterations, norm(bar_x - best_x), 'o')
%     hold on
    nexttile(3)
    title('|UB - LB|')
    plot(num_iterations, abs(UB-LB), 'o')
    hold on
    
end
hold off

if (num_iterations == max_iter)
    disp('Maximum iterations exceeded');
    exitFlag = 2;
    x_best = bar_x;
    return
end

end

