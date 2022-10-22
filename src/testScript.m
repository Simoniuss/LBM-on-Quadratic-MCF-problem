[Q, q, E, b, u] = loadMCF('MCF_100.mat');
[U, S, V] = compactSVD(E, size(E,1)-1);