function [mu, theta, l] = LBM(f, mu, l, B_z, B_alpha, lambda, bestl, m)
% Write some descriptions
%
%
%

if bestl
    
    [~, v, exitflag] = unstabilizedMP(B_z, B_alpha);

    % 0 means solution found, 2 means unbounded
    if exitflag == 0
        % Solution found, compute the best l
        l = lambda*f(mu) + (1- lambda)*v;
        [theta, exitflag] = thetaDP(B_z, B_alpha, l);
        if exitflag ~= 0
            error('Error in computing theta LBM bestl');
        end
    else
        % bestl unbounded
        if exitflag == 2
            emptyMP = true;
            while emptyMP % POTREBBE ESSERE NECESSARIO METTERE UN LIMITE SULLE ITERAZIONI
                [theta, exitflag] = thetaDP(B_z, B_alpha, l);
                
                % Solution found with the current l
                if exitflag == 0
                    emptyMP = false;
                % DP unbounded, change l
                elseif exitflag == 2
                    l = lambda*f(mu) - (1- lambda)*l;
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
            l = lambda*f(mu) - (1- lambda)*l;
        else
            error('Error in computing theta LBM arbirtrary l');
        end
    end
end

d = - B_z*theta;
v = max(B_z'*(mu + d) - B_alpha);

if f(mu+d) - f(mu) <= m*(v - f(mu))
    % SS (if not NS)
    disp('SS')
    disp(['d: ',num2str(norm(d))]);
    mu = mu+d;
else
    disp('NS');
end
      
end

