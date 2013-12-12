% Example script to generate synthetic data and run SPG method

% set sample size n, number of inputs p, number of outputs q
n = 100;
p = 30;
q = 15;

% generate synthetic data
[X,Y,Psi,Theta,B,G_x,G_y] = generate_synth_data(n,p,q);

% set lambda, gamma1 (for inputs), gamma2 (for outputs)
lambda = 2;
gamma1 = 1.5;
gamma2 = 1.5;

% run SPG for graph-structured input-output lasso with fused lasso penalty
BhatF = graph_siol_spg(X,Y,[],[],lambda,gamma1,gamma2,'fused');

% run SPG for graph-structured input-output lasso with group lasso penalty
BhatG = graph_siol_spg(X,Y,[],[],lambda,gamma1,gamma2,'group');

% plot results
GxTick = find(G_x > [0;G_x(1:end-1)]);
GyTick = find(G_y > [0;G_y(1:end-1)]);
minC = min(min(min(BhatF)),min(min(BhatG)));
maxC = max(max(max(BhatF)),max(max(BhatG)));
h = figure; hold on;
set(h,'Position',[10 10 1250 450]);
subplot('Position',[.02 .25 .06 .7]); imagesc(G_x); 
set(gca,'XTick',[],'YTick',GxTick-0.5,'YTickLabel',num2str(GxTick)); 
xlabel('G_X','FontSize',14);
subplot('Position',[.11 .25 .25 .7]); imagesc(B); 
set(gca,'XTick',GyTick-0.5,'XTickLabel',num2str(GyTick),'YTick',GxTick-0.5,'YTickLabel',num2str(GxTick));
title('True Regression Weights B','FontSize',14,'FontWeight','bold');
subplot('Position',[.39 .25 .25 .7]); imagesc(BhatF); caxis([minC,maxC]);
set(gca,'XTick',GxTick-0.5,'XTickLabel',num2str(GxTick),'YTick',GxTick-0.5,'YTickLabel',num2str(GxTick));
title('B Estimated with Fused Penalty','FontSize',14,'FontWeight','bold');
subplot('Position',[.67 .25 .25 .7]); imagesc(BhatG); caxis([minC,maxC]);
set(gca,'XTick',GyTick-0.5,'XTickLabel',num2str(GyTick),'YTick',GyTick-0.5,'YTickLabel',num2str(GyTick));
title('B Estimated with Group Penalty','FontSize',14,'FontWeight','bold');
colorbar('Position',[.95 .25 .03 .7]);
subplot('Position',[.11 .05 .25 .139]); imagesc(G_y'); 
set(gca,'XTick',GyTick-0.5,'XTickLabel',num2str(GyTick),'YTick',[]);
ylabel('G_Y','FontSize',14);

