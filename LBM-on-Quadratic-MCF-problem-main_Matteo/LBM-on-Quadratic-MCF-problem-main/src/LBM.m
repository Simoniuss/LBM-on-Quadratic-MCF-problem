function [mu, theta, l, d, steptype] = LBM(f, x, mu, l, B_z, B_alpha, lambda, bestl, m)
% Write some descriptions
%
%
%

if bestl
    
    [~, v, exitflag] = unstabilizedMP(B_z, B_alpha);

    % 0 means solution found, 2 means unbounded
    if exitflag == 0
        % Solution found, compute the best l
        l = lambda*f(x,mu) + (1 - lambda)*v;
        [theta, exitflag] = thetaDP(B_z, B_alpha, l);
        if exitflag ~= 0
            error('Error in computing theta LBM bestl');
        end
    else
        % bestl unbounded
        if exitflag == 2
            emptyMP = true;
            while emptyMP % POTREBBE ESSERE NECESSARIO METTERE UN LIMITE SULLE ITERAZIONI
                [theta, exitflag] = thetaDP(B_z, B_alpha, l, mu, f(x,mu));
                
                % Solution found with the current l
                if exitflag == 0
                    emptyMP = false;
                % DP unbounded, change l
                elseif exitflag == 2
                    l = lambda*f(x,mu) - (1 - lambda)*l;
                else
                    error('Error in computing theta LBM bestl unbounded');
                end
            end
        else
            error('Error in computing (mu,v) LBM bestl');
        end
    end
% arbitrary l
else
    emptyMP = true;
    while emptyMP
        [theta, exitflag] = thetaDP(B_z, B_alpha, l);
        % Solution found with the current l
        if exitflag == 0
            emptyMP = false;
            % DP unbounded, change l
        elseif exitflag == 2
            l = lambda*f(x,mu) - (1 - lambda)*l;
        else
            error('Error in computing theta LBM arbirtrary l');
        end
    end
end


d =  -B_z*theta;
v = max(B_z'*(mu + d) - B_alpha');
%d1 = findD(B_z, B_alpha, l);

%disp(norm(d))
%disp(norm(d1))

%disp(f(x, mu+d) - f(x,mu))
%disp(m*(v - f(x,mu)))


if f(x, mu+d) - f(x,mu) <= m*(v - f(x,mu))
    % SS (if not NS)
    steptype = 'SS';
    mu = mu+d;
else
    steptype = 'NS';
end
      
end

