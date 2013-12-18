function genC

subjects = 'ABCD';

p = 306;
q = 229;

% construct vectorized B indices
B_inds = zeros(p,q);
ind = 1;
for k = 1:q
    for j = 1:p
        B_inds(j,k) = ind;
        ind = ind + 1;
    end
end

C = cell(length(subjects), 1);
for s = 1:length(subjects)
    load(['./brainData/' subjects(s) '_sparse_data.mat']);
    C{s} = construct_fused_C_in(q, p, Psi, B_inds);
    fprintf('Subject %s complete\n', subjects(s));
    temp = sum(C{s}~=0,2);
    if min(temp)==2 && max(temp)==2
        fprintf('Current C is good to go!\n');
    else
        error('Rows of C do not have exactly 2 nonzero entries\n');
    end
end

save ./brainData/C_matrices_sparse_all.mat C;
end

% Function to generate matrix of size d1|E| x d1d2 that encodes a
% graph-guided fused lasso penalty
function C = construct_fused_C_in(q,p,Psi,B_inds)

numE = sum(sum(triu(Psi,1)~=0));
rowInd = zeros(q*numE,3);
colInd = zeros(q*p,2);
C = sparse(q*numE,q*p);
e = 1;
for l = 1:p
    for m = l+1:p
        if Psi(l,m) == 0; continue; end
        for k = 1:q
            C(numE*(k-1)+e,B_inds(l,k)) = abs(Psi(l,m));
            C(numE*(k-1)+e,B_inds(m,k)) = -(Psi(l,m));
            rowInd(numE*(k-1)+e,:) = [k l m];
            colInd(p*(k-1)+l,:) = [k l];
            colInd(p*(k-1)+m,:) = [k m];
        end
        e = e + 1;
    end
end

end
