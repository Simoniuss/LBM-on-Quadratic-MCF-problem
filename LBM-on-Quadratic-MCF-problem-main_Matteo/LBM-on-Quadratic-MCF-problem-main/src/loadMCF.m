function [Q, q, E, b, u] = loadMCF(matFilename)
% Take the mat file containing matrices and vectors of a generated MCF
% problem.
%       min 0.5 * x' * Q * x + q * x
%       s.t. Ex=b, 0 <= x <= u
% Q: diagonal of a matrix postive semi-definite Nx1
% q: vector Nx1
% E: incidence matrix of a connected digraph MxN
% b: vector Mx1 for Ex=b
% u: vector for right box constraints Nx1
MCF = load(matFilename);
Q = MCF.Q';
q = MCF.q';
E = full(MCF.E);
b = MCF.b';
u = MCF.u';
end

