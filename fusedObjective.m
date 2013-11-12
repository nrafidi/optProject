function obj = fusedObjective(X, Y, B, slackVar, lambda, gamma1, gamma2, Psi, Theta)

p = size(Psi, 1);
q = size(Theta, 1);

obj = trace((Y - X*B)*(Y - X*B)') + lambda*sum(sum((B.^2)./slackVar{1}));

if gamma1 ~= 0
    pen1 = 0;
    for l = 1:p
        for m = (l+1):p
            pen1 = pen1 + Psi(l,m)^2*sum((B(l,:) - sign(Psi(l,m))*B(m,:)).^2./squeeze(slackVar{2}(:, l, m)'));
        end
    end
end
if gamma2 ~= 0
    pen2 = 0;
    for r = 1:q
        for s = (r+1):q
            pen2 = pen2 + Theta(r,s)^2*sum((B(:,r) - sign(Theta(r,s))*B(:,s)).^2./squeeze(slackVar{3}(:, r, s)));
        end
    end
end
if gamma1 ~= 0
    obj = obj + gamma1*pen1;
end
if gamma2 ~= 0
    obj = obj + gamma2*pen2;
end

end