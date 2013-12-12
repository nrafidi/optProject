clear;

K=10; % dimension of Y
J=10;  % dimension of X 
N=100; % sample size

option.cortype=2;
option.nE=5*J;  % pre-specify no of edges for constructing the graph, (another method, threshold
%correlation, has also been implemented in gennetwork.m

option.maxiter=10000;
option.tol=1e-8;  
option.verbose=true;
option.display_iter=50;
option.mu=1e-4;  % smaller mu leads more accurate solution


gamma=10;
lambda=gamma;

[Y, X, B]=simulate_data(N,J,K);
[C, CNorm, E, Ecoef, Esign, R]=gennetwork(Y, option);

[beta, obj, density, iter, time]=SPG_multi(Y, X, gamma, lambda, C, CNorm, option);
%[beta, obj, density, iter,  time]=SPG_multi_linesearch(Y, X, gamma, lambda, C,  option);