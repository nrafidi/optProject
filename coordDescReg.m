function [B, slackVar, lambda, gamma1, gamma2] = coordDescReg(X, Y, coordFuncs, objFunc, Psi, Theta)

%Parameter Selection - rough iterative grid search
regParams = [.0000001 .000001 .00001 .0001 .001 .01 .1 .5 1 5 10 50 100 500 1000 10000 20000 50000 ...
    100000 500000 1000000 5000000 10000000];
R = length(regParams);

n = size(X,1);
k = n; %number of folds for CV
cvInd = crossvalind('Kfold', n, k);

%Choose Lambda, fixing gamma1 and gamma2 to 0
cverrs = zeros(R, 1);
for r = 1:R
    %Cross validation
    for kk = 1:k
        Breg = coordDesc(X(cvInd~=kk, :), Y(cvInd~=kk,:), regParams(r), 0, 0, coordFuncs, objFunc, Psi, Theta);
        cverrs(r) = cverrs(r) + norm(Y(cvInd==kk,:) - X(cvInd==kk,:)*Breg)/k;
    end
end
[~, ind] = min(cverrs);
lambda = regParams(ind);

%Choose gamma1
cverrs = zeros(r, 1);
for r = 1:R
    %Cross validation
    for kk = 1:k
        Breg = coordDesc(X(cvInd~=kk, :), Y(cvInd~=kk,:), lambda, regParams(r), 0, coordFuncs, objFunc, Psi, Theta);
        cverrs(r) = cverrs(r) + norm(Y(cvInd==kk,:) - X(cvInd==kk,:)*Breg)/k;
    end
end
[~, ind] = min(cverrs);
gamma1 = regParams(ind);

%Choose gamma2
cverrs = zeros(r, 1);
for r = 1:R
    %Cross validation
    for kk = 1:k
        Breg = coordDesc(X(cvInd~=kk, :), Y(cvInd~=kk,:), lambda, gamma1, regParams(r), coordFuncs, objFunc, Psi, Theta);
        cverrs(r) = cverrs(r) + norm(Y(cvInd==kk,:) - X(cvInd==kk,:)*Breg)/k;
    end
end
[~, ind] = min(cverrs);
gamma2 = regParams(ind);

%Parameter Selection - gradient descent -
tolerance = 1e-4;
step = 1e-2;
h = 0.1; %I have legit no idea what to make this

[B, slackVar] = coordDesc(X, Y, lambda, gamma1, gamma2, coordFuncs, objFunc, Psi, Theta);
currObj = feval(objFunc, X, Y, B, slackVar, lambda, gamma1, gamma2, Psi, Theta);
prevObj = 1e10;
while abs(currObj - prevObj) > tolerance
    %Estimate gradient:
    [lambdaB, lamSlack] = coordDesc(X, Y, lambda + h, gamma1, gamma2, coordFuncs, objFunc, Psi, Theta);
    [gamma1B, gam1Slack] = coordDesc(X, Y, lambda, gamma1+h, gamma2, coordFuncs, objFunc, Psi, Theta);
    [gamma2B, gam2Slack] = coordDesc(X, Y, lambda, gamma1, gamma2+h, coordFuncs, objFunc, Psi, Theta);
    
    delLambda = - step*(feval(objFunc, X, Y, lambdaB, lamSlack, lambda, gamma1, gamma2, Psi, Theta) - currObj)/h;
    delGamma1 = - step*(feval(objFunc, X, Y, gamma1B, gam1Slack, lambda, gamma1, gamma2, Psi, Theta) - currObj)/h;
    delGamma2 = - step*(feval(objFunc, X, Y, gamma2B, gam2Slack, gamma1, gamma2, Psi, Theta) - currObj)/h;
    
    lambda = lambda + delLambda;
    gamma1 = gamma1 + delGamma1;
    gamma2 = gamma2 + delGamma2;
    
    prevObj = currObj;
    [B, slackVar] = coordDesc(X, Y, lambda, gamma1, gamma2, coordFuncs, objFunc, Psi, Theta);
    currObj = feval(objFunc, X, Y, B, slackVar, lambda, gamma1, gamma2, Psi, Theta);
end

%The last B learned in the above procedure is the one we return.

end