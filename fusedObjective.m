function obj = fusedObjective(X, Y, B, slackVar, lambda, gamma1, gamma2, Psi, Theta)

p = size(Psi, 1);
q = size(Theta, 1);

obj = trace((Y - X*B)'*(Y - X*B)) + lambda*sum(sum((B.^2)./slackVar{1}));

pen1 = 0;
for l = 1:p
    for m = l:p
        pen1 = pen1 + Psi(l,m)^2*sum((B(l,:) - sign(Psi(l,m))*B(m,:))./squeeze(slackVar{2}(:, l, m)'));
    end
end

pen2 = 0;
for r = 1:q
    for s = r:q
        pen2 = pen2 + Theta(r,s)^2*sum((B(:,r) - sign(Theta(r,s))*B(:,s))./squeeze(slackVar{3}(:, r, s)));
    end
end

obj = obj + gamma1*pen1 + gamma2*pen2;

end