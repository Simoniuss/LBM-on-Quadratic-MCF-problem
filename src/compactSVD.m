function [U, S, V] = compactSVD(E, r)
% Return the compact SVD of the matrix E which is the matrix truncated to
% the rank of the matrix r. If the rank is not specified, it is computed 
% by the function.
% E: MxN
% r: rank of E
% U: Mxr
% S: rxr
% V: NxR
if ~exist('r', 'var')
    r = rank(E);
end
[U, S, V] = svd(E, 'econ');
U = U(:, 1:r);
S = S(1:r, 1:r);
V = V(:, 1:r);
end

