function [X,Y,Psi,Theta,B,G_x,G_y] = generate_synth_data(n,p,q)

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
Psi(abs(Psi) < 0.2) = 0;

% calculate correlation matrix for the outputs
Theta = Y'*Y/(n-1);
Theta(abs(Theta) < 0.2) = 0;

% plot results
GxTick = find(G_x > [0;G_x(1:end-1)]);
GyTick = find(G_y > [0;G_y(1:end-1)]);
minC = min(min(min(Psi)),min(min(Theta)));
maxC = max(max(max(Psi)),max(max(Theta)));
h = figure; hold on;
set(h,'Position',[10 10 1250 450]);
subplot('Position',[.02 .25 .06 .7]); imagesc(G_x); 
set(gca,'XTick',[],'YTick',GxTick-0.5,'YTickLabel',num2str(GxTick)); 
xlabel('G_X','FontSize',14);
subplot('Position',[.11 .25 .25 .7]); imagesc(B); 
set(gca,'XTick',GyTick-0.5,'XTickLabel',num2str(GyTick),'YTick',GxTick-0.5,'YTickLabel',num2str(GxTick));
title('True Regression Weights B','FontSize',14,'FontWeight','bold');
subplot('Position',[.39 .25 .25 .7]); imagesc(Psi); caxis([minC,maxC]);
set(gca,'XTick',GxTick-0.5,'XTickLabel',num2str(GxTick),'YTick',GxTick-0.5,'YTickLabel',num2str(GxTick));
title('Input Structure \Psi','FontSize',14);
subplot('Position',[.67 .25 .25 .7]); imagesc(Theta); caxis([minC,maxC]);
set(gca,'XTick',GyTick-0.5,'XTickLabel',num2str(GyTick),'YTick',GyTick-0.5,'YTickLabel',num2str(GyTick));
title('Output Structure \Theta','FontSize',14);
colorbar('Position',[.95 .25 .03 .7]);
subplot('Position',[.11 .05 .25 .139]); imagesc(G_y'); 
set(gca,'XTick',GyTick-0.5,'XTickLabel',num2str(GyTick),'YTick',[]);
ylabel('G_Y','FontSize',14);
% subplot('Position',[.39 .05 .25 .139]); imagesc(G_y'); 
% set(gca,'XTick',GyTick-0.5,'XTickLabel',num2str(GyTick),'YTick',[]);
% subplot('Position',[.67 .05 .25 .139]); imagesc(G_y'); 
% set(gca,'XTick',GyTick-0.5,'XTickLabel',num2str(GyTick),'YTick',[]);

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
    for c = connInd
        while isempty(G_in==c), c = randsample(numGin,1); end
    end
    for c = connInd
        B(G_in==c,G_out==gout) = 0.8;
    end
end

end

