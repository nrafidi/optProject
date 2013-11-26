function [B, slackVar, lambda, gamma1, gamma2, totalIT, avgIT] = coordDescReg(X, Y, coordFuncs, objFunc, Psi, Theta)

%Parameter Selection - rough iterative grid search
lambdas = fliplr([1e-5 1e-4 0.001, 0.01, 0.1, 1, 10, 100]);
gammas = fliplr([1e-5 1e-4 0.001, 0.01, 0.1, 1, 10, 100]);
R = length(gammas);

n = size(X,1);
k = 2; %number of folds for CV
cvInd = crossvalind('Kfold', n, k);

% lambda = 0;
gamma1 = 0;
gamma2 = 0;
p = size(X,2);
q = size(Y,2);
avgIT = 0;

%Choose Lambda, fixing gamma1 and gamma2 to 0
Breg = rand(p,q);
cverrs = zeros(R, 1);
numIT = zeros(R, k);
for r = 1:R
    %Cross validation
    for kk = 1:k
        [Breg, ~, numIT(r,kk)] = coordDesc(X(cvInd~=kk, :), Y(cvInd~=kk,:), lambdas(r), gamma1, gamma2, coordFuncs, objFunc, Psi, Theta, Breg);
        cverrs(r) = cverrs(r) + norm(Y(cvInd==kk,:) - X(cvInd==kk,:)*Breg)/k;
    end
%     fprintf('Woof\n');
end
[~, ind] = min(cverrs);
avgIT = avgIT + mean(numIT(ind,:),2);
lambda = lambdas(ind);

% Plotting Code
% B = coordDesc(X, Y, lambda, gamma1, gamma2, coordFuncs, objFunc, Psi, Theta, rand(p,q));
% subplot(2,2,2);
% imagesc(B)
% colorbar;
% title(sprintf('\\lambda = %d, \\gamma_1 = %d, \\gamma_2 = %d', lambda, gamma1, gamma2));
fprintf('Lambda done\n');

%Choose gamma1
Breg = rand(p,q);
cverrs = zeros(R, 1);
numIT = zeros(R, k);
for r = 1:R
    %Cross validation
    for kk = 1:k
        [Breg, ~, numIT(r,kk)] = coordDesc(X(cvInd~=kk, :), Y(cvInd~=kk,:), lambda, gammas(r), gamma2, coordFuncs, objFunc, Psi, Theta, Breg);
        cverrs(r) = cverrs(r) + norm(Y(cvInd==kk,:) - X(cvInd==kk,:)*Breg)/k;
    end
%     fprintf('Meow\n');
end
[~, ind] = min(cverrs);
avgIT = avgIT + mean(numIT(ind,:),2);
gamma1 = gammas(ind);

%Plotting code
% B = coordDesc(X, Y, lambda, gamma1, gamma2, coordFuncs, objFunc, Psi, Theta, rand(p,q));
% subplot(2,2,3);
% imagesc(B)
% colorbar;
% title(sprintf('\\lambda = %d, \\gamma_1 = %d, \\gamma_2 = %d', lambda, gamma1, gamma2));
fprintf('gamma1 done\n');

%Choose gamma2
Breg = rand(p,q);
cverrs = zeros(R, 1);
numIT = zeros(R, k);
for r = 1:R
    %Cross validation
    for kk = 1:k
        [Breg, ~, numIT(r,kk)] = coordDesc(X(cvInd~=kk, :), Y(cvInd~=kk,:), lambda, gamma1, gammas(r), coordFuncs, objFunc, Psi, Theta, Breg);
        cverrs(r) = cverrs(r) + norm(Y(cvInd==kk,:) - X(cvInd==kk,:)*Breg)/k;
    end
%     fprintf('Oink\n');
end
[~, ind] = min(cverrs);
avgIT = avgIT + mean(numIT(ind,:),2);
gamma2 = gammas(ind);
[B, slackVar, totalIT] = coordDesc(X, Y, lambda, gamma1, 0, coordFuncs, objFunc, Psi, Theta, rand(p,q));

%Plotting Code
% subplot(2,2,4);
% imagesc(B)
% colorbar;
% title(sprintf('\\lambda = %d, \\gamma_1 = %d, \\gamma_2 = %d', lambda, gamma1, gamma2));
% suptitle('B after Different Stages of Initial Parameter Search');
% keyboard
fprintf('gamma2 done\n');

% Parameter Selection - gradient descent -
% tolerance = 1e-4;
% step = 1e-6;
% h = 0.1; %I have legit no idea what to make this
% 
% lambdaB = rand(p,q);
% gamma1B = lambdaB;
% gamma2B = lambdaB;
% [B, slackVar] = coordDesc(X, Y, lambda, gamma1, gamma2, coordFuncs, objFunc, Psi, Theta, rand(p,q));
% currObj = feval(objFunc, X, Y, B, slackVar, lambda, gamma1, gamma2, Psi, Theta);
% prevObj = 1e10;
% while abs(currObj - prevObj) > tolerance
%     %Estimate gradient:
%     [lambdaB, lamSlack] = coordDesc(X, Y, lambda + h, gamma1, gamma2, coordFuncs, objFunc, Psi, Theta, lambdaB);
%     [gamma1B, gam1Slack] = coordDesc(X, Y, lambda, gamma1+h, gamma2, coordFuncs, objFunc, Psi, Theta, gamma1B);
%     [gamma2B, gam2Slack] = coordDesc(X, Y, lambda, gamma1, gamma2+h, coordFuncs, objFunc, Psi, Theta, gamma2B);
%     
%     delLambda = - step*(feval(objFunc, X, Y, lambdaB, lamSlack, lambda, gamma1, gamma2, Psi, Theta) - currObj)/h;
%     delGamma1 = - step*(feval(objFunc, X, Y, gamma1B, gam1Slack, lambda, gamma1, gamma2, Psi, Theta) - currObj)/h;
%     delGamma2 = - step*(feval(objFunc, X, Y, gamma2B, gam2Slack, lambda, gamma1, gamma2, Psi, Theta) - currObj)/h;
%     
%     lambda = lambda + delLambda;
%     gamma1 = gamma1 + delGamma1;
%     gamma2 = gamma2 + delGamma2;
%     
%     prevObj = currObj;
%     [B, slackVar] = coordDesc(X, Y, lambda, gamma1, gamma2, coordFuncs, objFunc, Psi, Theta, B);
%     currObj = feval(objFunc, X, Y, B, slackVar, lambda, gamma1, gamma2, Psi, Theta);
% end
% 
% %The last B learned in the above procedure is the one we return.
% figure
% imagesc(B);
% title(sprintf('\\lambda = %d, \\gamma_1 = %d, \\gamma_2 = %d', lambda, gamma1, gamma2));
end