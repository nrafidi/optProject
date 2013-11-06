function B = fusedFusedUpdate(X, Y, B, Psi, Theta, slackVars, lambda, gamma1, gamma2)

p = size(X,2);
q = size(Y,2);
indX = 1:p;

% newB = zeros(p,q); tralalalala nothing is atomic

for j = 1:p
    for k = 1:q
        B(j,k) = sum(X(:,j).*(Y(:,k) - X(:,indX ~=j)*B(indX ~=j, k)));
        B(j,k) = B(j,k) + gamma1*sum((Psi(j,:).^2.*sign(Psi(j,:)).*B(:,k))./squeeze(slackVars{2}(j, k, :)));
        B(j,k) = B(j,k) + gamma2*sum((Theta(k,:).^2.*sign(Theta(k,:)).*B(j,:))./squeeze(slackVars{3}(j, k, :)));
        B(j,k) = B(j,k)/(sum(X(:,j).^2) + lambda/slackVars{1}(j,k) + ...
            gamma1*sum(Psi(j,:).^2./squeeze(slackVars{2}(j, k, :))) + ...
            gamma2*sum(Theta(k,:).^2./squeeze(slackVars{3}(j, k, :))));
    end
end

end