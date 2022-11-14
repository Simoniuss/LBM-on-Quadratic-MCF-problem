[Q, q, E, b, u] = loadMCF('MCF_1000.mat');
%[U, S, V, um] = compactSVD(E);
%simpleProblem;

epsilon = 1e-4;
%l = 6.3282e+03;
l=3.4687e+08;
lambda = 0.5;
best_l = false;
m = 0.1;
max_iter = 50;

[x_best, exitFlag] = QMCF_solver(Q, q, E, b, u, epsilon, l, lambda, best_l, m, max_iter);