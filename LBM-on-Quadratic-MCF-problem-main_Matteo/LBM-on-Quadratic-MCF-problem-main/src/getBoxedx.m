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

