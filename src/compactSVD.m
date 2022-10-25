function [U, S, V, um] = compactSVD(E)
% Return the compact SVD of the matrix E, which is the incidence matrix of
% a weakly connected component graph, truncated to the rank of the matrix 
% r=m-1. It returns also the cokernel(E) which is U_m column.
%
% E: MxN
% U: Mx(M-1)
% S: (M-1)x(M-1)
% V: Nx(M-1)
% um: Mx1
m = size(E,1);
[U, S, V] = svd(E, 'econ');
um = U(:,end);
r = m-1;
U = U(:, 1:r);
S = S(1:r, 1:r);
V = V(:, 1:r);
end

