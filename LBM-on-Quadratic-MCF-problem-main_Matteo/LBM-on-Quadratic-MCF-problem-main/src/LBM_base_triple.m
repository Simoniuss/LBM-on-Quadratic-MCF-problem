function [d, exitflag, theta, f_B] = LBM_base_triple(B_z, B_f, B_mu, l, mu);

muSize = size(mu,1);

MUPROB = false;

if(MUPROB)

    %MU PROB
    muVar = sdpvar(muSize,1);
    %try to rewrite it
    lin = [];
    for i = 1:size(B_z,2)
        lin = [lin ( B_f(i)+ B_z(:,i)'*(muVar- B_mu(:,i)) )];
    end
    
    Constraints = [l  >= max(lin)];
    
    %MU PRob
    Objective = (norm(muVar - mu)^2) /2;
    
    %options = sdpsettings('verbose', 0, 'solver','quadprog');
    options = sdpsettings('verbose', 0, 'solver','mosek');
    %options = sdpsettings('verbose', 0);
    
    sol = optimize(Constraints, Objective, options);
    
    %MU PROB
    mu = value(muVar);
    for i = 1:size(B_z,2)
        if(i==1)
            f_B = ( B_f(i)+ B_z(:,i)'*(mu- B_mu(:,i)) );
        else
            other = ( B_f(i)+ B_z(:,i)'*(mu- B_mu(:,i)) );
            if(other > f_B)
                f_B = other;
            end
        end
    end

else

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %D prob
    dVar = sdpvar(muSize,1);
    
    %try to rewrite it
    lin = [];
    for i = 1:size(B_z,2)
        lin = [lin ( B_f(i)+ B_z(:,i)'*((mu+dVar)- B_mu(:,i)) )];
    end
    
    %disp(f_B)
    Constraints = [l  >= max(lin)];
    
    %D PROB
    Objective = (norm(dVar)^2) / 2;
    
    %options = sdpsettings('verbose', 0, 'solver','quadprog');
    options = sdpsettings('verbose', 0, 'solver','mosek');
    %options = sdpsettings('verbose', 0);
    
    sol = optimize(Constraints, Objective, options);
    
    %D PROB
    d = value(dVar);
    for i = 1:size(B_z,2)
        if(i==1)
            f_B = ( B_f(i)+ B_z(:,i)'*((mu+d)- B_mu(:,i)) );
        else
            other = ( B_f(i)+ B_z(:,i)'*((mu+d)- B_mu(:,i)) );
            if(other > f_B)
                f_B = other;
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%

theta = dual(Constraints);
exitflag = sol.problem;

end
