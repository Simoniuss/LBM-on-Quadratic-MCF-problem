[Q, q, E, b, u] = loadMCF('MCF_1000.mat');
[U, S, V, um] = compactSVD(E);

epsilon = 1.0e-6;
l = -0.01;
lambda = 0.5;
best_l = false;
m = 0.5;
max_iter = 100;

[x_best, exitFlag] = QMCF_solver(Q, q, E, b, u, epsilon, l, lambda, best_l, m, max_iter);