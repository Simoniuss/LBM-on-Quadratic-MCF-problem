%[Q, q, E, b, u] = loadMCF('MCF_1000.mat');
%[U, S, V, um] = compactSVD(E);
simpleProblem;

u(1) = 1;

epsilon = 1e-03;
%l = 6.3282e+03;
l = -6.3282e+02;
%l = 3.4687e+08;
lambda = 0.005;
%lambda = 0.05;
best_l = false;
m = 0.1;
max_iter = 5000;

[x_best, exitFlag] = QMCF_solver_v4(Q, q, E, b, u, epsilon, l, lambda, best_l, m, max_iter);