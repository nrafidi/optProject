function runBrainExpSPG
subjects = 'ABCD';
load ../brainData/semantic_data.mat;
load ../brainData/C_matrices_sparse_all.mat;
if exist('../brainData/C_out_all.mat', 'file');
load ../brainData/C_out_all.mat;
else
    p = 306;
    q = 218;
    % construct vectorized B indices
    B_inds = zeros(p,q);
    ind = 1;
    for k = 1:q
        for j = 1:p
            B_inds(j,k) = ind;
            ind = ind + 1;
        end
    end
    C_out = construct_fused_C_out(p,q,Theta(1:218,1:218),B_inds);
    temp = sum(C_out~=0,2);
    if min(temp)==2 && max(temp)==2
        fprintf('Current C is good to go!\n');
    else
        error('Rows of C do not have exactly 2 nonzero entries\n');
    end
save ../brainData/C_out_all.mat C_out_all;
end
addpath ./spg_package/SPG_singletask/
cd ./spg_package/SPG_singletask/
dir
install_mex;
cd ../../

%% SPG Code - Training
penFlag = 'fused';


betas = cell(length(subjects), 1);
for s = 1:length(subjects)
    load(['../brainData/' subjects(s) '_sparse_data.mat']);
    betas{s} = graph_siol_spg_reg(X, Y(:, 1:218), C{s}, C_out, penFlag);
    fprintf('Subject %s complete\n', subjects(s));
end

save ../brainData/betasSPG_sparse_all.mat betas;
end


% Function to generate matrix of size d1|E| x d1d2 that encodes a
% graph-guided fused lasso penalty
function C = construct_fused_C_out(p,q,Theta,B_inds)

numE = sum(sum(triu(Theta,1)~=0));
rowInd = zeros(p*numE,3);
colInd = zeros(p*q,2);
C = sparse(p*numE,p*q);
e = 1;
for r = 1:q
    for s = r+1:q
        if Theta(r,s) == 0; continue; end
        for j = 1:p
            C(numE*(j-1)+e,B_inds(j,r)) = abs(Theta(r,s));
            C(numE*(j-1)+e,B_inds(j,s)) = -(Theta(r,s));
            rowInd(numE*(j-1)+e,:) = [j r s];
            colInd(q*(j-1)+r,:) = [j r];
            colInd(q*(j-1)+s,:) = [j s];
        end
        e = e + 1;
    end
end

end
