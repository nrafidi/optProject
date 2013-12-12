% Function to run smoothing proximal gradient descent for graph-guided 
% structured input-output lasso, using either the fused lasso penalty or 
% the group lasso penalty

% Inputs:
%   X: n-by-p matrix of input features, one feature per column
%   Y: n-by-q matrix of output tasks, one task per column
%   Psi: p-by-p matrix encoding graph structure among inputs [optional]
%   Theta: q-by-q matrix encoding graph structure among outputs [optional]
%   lambda: weight for L1 norm penalty [optional]
%   gamma1: weight for input structured sparsity penalty [optional]
%   gamma2: weight for output structured sparsity penalty [optional]
%   penFlag: indicates whether to use fused lasso penalty ('fused') or 
%       group lasso penalty ('group') [optional]

function Bhat = graph_siol_spg(X,Y,Psi,Theta,lambda,gamma1,gamma2,penFlag)

% get dimensions
[n,p] = size(X);
[~,q] = size(Y);

% standardize X and Y
X = zscore(X);
Y = zscore(Y);

% compute Psi and Theta if not provided
if ~exist('Psi','var') || isempty(Psi)
    Psi = X'*X/(n-1);
    Psi(abs(Psi) < 0.2) = 0;
end
if ~exist('Theta','var') || isempty(Theta)
    Theta = Y'*Y/(n-1);
    Theta(abs(Theta) < 0.2) = 0;
end

% set lambda and gammas if not provided
if ~exist('lambda','var') || isempty(lambda)
    lambda = 1;
end
if ~exist('gamma1','var') || isempty(gamma1)
    gamma1 = 1;
end
if ~exist('gamma2','var') || isempty(gamma2)
    gamma2 = 1;
end

% set penalty flag if not provided
if ~exist('penFlag','var') || isempty(penFlag)
    penFlag = 'fused';
elseif ~strcmp(penFlag,'fused') && ~strcmp(penFlag,'group')
    error('Specify fused or group penalty\n');
end

% set options for optimization
option.maxiter = 10000;
option.tol = 1e-4;  
option.verbose = true;
option.display_iter = 50;
option.mu = 1e-4;

% construct C matrix
if strcmp(penFlag,'fused')
    C_in = construct_fused_C(q,p,Psi);
    C_out = construct_fused_C(p,q,Theta);
    C = [gamma1*C_in ; gamma2*C_out];
elseif strcmp(penFlag,'group')
    [C_in,G_in] = construct_group_C(q,p,Psi);
    [C_out,G_out] = construct_group_C(p,q,Theta);
    C = [gamma1*C_in ; gamma2*C_out];
    G = [G_in ; G_out];
end 

% reshape data matrices X and Y
y = Y(:);
x = sparse(n*q,p*q);
for k = 1:q
    starti = (k-1)*n;
    startj = (k-1)*p;
    x(starti+1:starti+n,startj+1:startj+p) = X;
end

% run smoothing proximal gradient descent
if strcmp(penFlag,'fused')    
    fprintf('Running SPG method for optimizing SIOL-Fused objective...\n');
    [beta,obj,density,iter,time] = SPG_linesearch('graph',y,x,1,lambda,C,option);
elseif strcmp(penFlag,'group')
    fprintf('Running SPG method for optimizing SIOL-Group objective...\n');
    [beta,obj,density,iter,time] = SPG_linesearch('group',y,x,1,lambda,C,option,G);
end

% reshape estimated parameters beta
Bhat = reshape(beta,p,q);

end

% Function to generate matrix of size d1|E| x d1d2 that encodes a
% graph-guided fused lasso penalty
function C = construct_fused_C(d1,d2,Sigma)


numE = sum(sum(triu(Sigma,1)~=0));
C = sparse(d1*numE,d1*d2);
e = 1;
for r = 1:d2
    for s = r+1:d2
        if Sigma(r,s) == 0; continue; end
        for j = 1:d1
            C(numE*(j-1)+e,d2*(j-1)+r) = abs(Sigma(r,s));
            C(numE*(j-1)+e,d2*(j-1)+s) = -(Sigma(r,s));
        end
        e = e + 1;
    end
end

end

% Function to generate matrix of size d1|E| x d1d2 that encodes a 
% graph-guided group lasso penalty
function [C,G] = construct_group_C(d1,d2,Sigma)

numE = sum(sum(triu(Sigma,1)~=0));
C = sparse(2*d1*numE,d1*d2);
G = zeros(d1*numE,3);
ind = 1;
e = 1;
for r = 1:d2
    for s = r+1:d2
        if Sigma(r,s) == 0; continue; end
        for j = 1:d1
            C(ind,d2*(j-1)+r) = abs(Sigma(r,s));
            C(ind+1,d2*(j-1)+s) = abs(Sigma(r,s));
            G(numE*(j-1)+e,:) = [ind ind+1 2];
            ind = ind + 2;
        end
        e = e + 1;
    end
end

end
