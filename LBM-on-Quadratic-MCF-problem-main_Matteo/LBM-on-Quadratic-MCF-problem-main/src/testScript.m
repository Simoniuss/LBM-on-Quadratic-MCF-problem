%[Q, q, E, b, u] = loadMCF('MCF_50.mat');
%[U, S, V, um] = compactSVD(E);
simpleProblem;

u(1) = 1;

epsilon = 1e-03;
%l = 6.3282e+03;
%l = -8.9282e+03;
l = -4.4687e+03;
lambda = 0.4;
%lambda = 0.05;
best_l = false;
m = 0.001;
max_iter = 10000;

[x_best, exitFlag] = QMCF_solver_Prof_mixed(Q, q, E, b, u, epsilon, l, lambda, best_l, m, max_iter);