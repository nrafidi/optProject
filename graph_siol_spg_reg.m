function Bhat = graph_siol_spg_reg(X,Y,C_in,C_out, penFlag)
%graph_siol_spg_reg performs 2-fold cv to choose hyperparameters for
%grap_siol_spg

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

%Choose Lambda
cverrs = zeros(R, 1);
for r = 1:R
    %Cross validation
    for kk = 1:k
        Breg = graph_siol_spg(X(cvInd~=kk, :), Y(cvInd~=kk,:), C_in, C_out, lambdas(r), gamma1, gamma2, penFlag);
        cverrs(r) = cverrs(r) + norm(Y(cvInd==kk,:) - X(cvInd==kk,:)*Breg)/k;
    end
end
[~, ind] = min(cverrs);
lambda = lambdas(ind);

%Choose Gamma1
cverrs = zeros(R, 1);
for r = 1:R
    %Cross validation
    for kk = 1:k
        Breg = graph_siol_spg(X(cvInd~=kk, :), Y(cvInd~=kk,:), C_in, C_out, lambda, gammas(r), gamma2, penFlag);
        cverrs(r) = cverrs(r) + norm(Y(cvInd==kk,:) - X(cvInd==kk,:)*Breg)/k;
    end
end
[~, ind] = min(cverrs);
gamma1 = gammas(ind);


%Choose Gamma2
cverrs = zeros(R, 1);
for r = 1:R
    %Cross validation
    for kk = 1:k
        Breg = graph_siol_spg(X(cvInd~=kk, :), Y(cvInd~=kk,:), C_in, C_out, lambda, gamma1, gammas(r), penFlag);
        cverrs(r) = cverrs(r) + norm(Y(cvInd==kk,:) - X(cvInd==kk,:)*Breg)/k;
    end
end
[~, ind] = min(cverrs);
gamma2 = gammas(ind);

Bhat = graph_siol_spg(X, Y, C_in, C_out, lambda, gamma1, gamma2, penFlag);
end

