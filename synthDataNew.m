function [X, Y, Psi, Theta, B] = synthDataNew(n, p, q)

% generate groups among the inputs
G_x = generateGroups(p);

% generate groups among the outputs
G_y = generateGroups(q);

% generate block sparse covariance matrix for the inputs
S_x = generateBlockCovariance(G_x);

% generate conditional covariance matrix for the outputs
S_y = eye(q);

% generate structured sparse regression coefficients
B = generateRegressionWeights(G_x,G_y);

% generate inputs X
X = mvnrnd(zeros(n,p),S_x);
X = zscore(X);

% generate outputs Y
Y = mvnrnd(X*B,S_y);
Y = zscore(Y);

% calculate correlation matrix for the inputs
Psi = X'*X/(n-1);

% calculate correlation matrix for the outputs
Theta = Y'*Y/(n-1);

end

% partition indices 1:v into groups
function G = generateGroups(numVar)

% fix number of groups
numGroup = floor(numVar/3);

% assign each variable to a group
groups = randi(numGroup,numVar,1);

% order variables according to groups
G = sort(groups);

end

% generate block sparse covariance matrix
function S = generateBlockCovariance(G)

% initialize S
v = length(G);
S = zeros(v,v);

% get number of groups
numG = length(G);

% construct each block of S
for g = 1:numG
    gInd = find(G==g);
    gLen = length(gInd);
    S(gInd,gInd) = ones(gLen,gLen);
end

end

% generate regression weights
function B = generateRegressionWeights(G_in,G_out)

% initialize B
p = length(G_in);
q = length(G_out);
B = zeros(p,q);

% get number of groups
numGin = max(G_in);
numGout = max(G_out);

% set distribution over number of connections
connDist = sort(4.^(0:numGin-1),'descend');
connDist = connDist./sum(connDist);

% assign each output group to one or more input groups
for gout = 1:numGout
    numConn = find(mnrnd(1,connDist));
    connInd = randsample(numGin,numConn)';
    display(num2str(numConn));
    for c = connInd
        B(G_in==c,G_out==gout) = 0.8;
    end
end

end

