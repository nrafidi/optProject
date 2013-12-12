function R=shrink(A, g_idx)
    V=size(g_idx,1);
    R=zeros(size(A));
    for v=1:V
        idx=g_idx(v,1):g_idx(v,2);
        num=g_idx(v,3);
        gnorm=sqrt(sum(A(idx,:).^2));
        gnorm(gnorm<1)=1;
        R(idx,:)=A(idx,:)./repmat(gnorm,  num, 1);        
    end
end