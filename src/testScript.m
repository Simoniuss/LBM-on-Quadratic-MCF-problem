[Q, q, E, b, u] = loadMCF('MCF_1000.mat');
%[U, S, V, um] = compactSVD(E);

epsilon = 1e-6;
l = -1e-10;
lambda = 0.5;
best_l = true;
m = 0.5;
max_iter = 100;

[x_best, exitFlag] = QMCF_solver(Q, q, E, b, u, epsilon, l, lambda, best_l, m, max_iter);