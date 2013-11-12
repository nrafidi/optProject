function B = fusedFusedUpdate(X, Y, B, Psi, Theta, slackVars, lambda, gamma1, gamma2)

p = size(X,2);
q = size(Y,2);
indX = 1:p;

% newB = zeros(p,q); tralalalala nothing is atomic

for j = 1:p
    for k = 1:q
        if isnan(B(j,k)) || isinf(B(j,k))
            keyboard;
        end
        Bold = B;
        B(j,k) = sum(X(:,j).*(Y(:,k) - X(:,indX ~=j)*Bold(indX ~=j, k)));
        if gamma1 ~= 0
            B(j,k) = B(j,k) + gamma1*sum((Psi(j,:)'.^2.*sign(Psi(j,:)').*Bold(:,k))./squeeze(slackVars{2}(k, j, :)));
        end
        if gamma2 ~= 0
            B(j,k) = B(j,k) + gamma2*sum((Theta(k,:)'.^2.*sign(Theta(k,:)').*Bold(j,:)')./squeeze(slackVars{3}(j, k, :)));
        end
        B(j,k) = B(j,k)/(sum(X(:,j).^2) + lambda/slackVars{1}(j,k) + ...
            gamma1*sum(Psi(j,:)'.^2./squeeze(slackVars{2}(k, j, :))) + ...
            gamma2*sum(Theta(k,:)'.^2./squeeze(slackVars{3}(j, k, :))));
        if isnan(B(j,k)) || isinf(B(j,k))
            keyboard;
        end
    end
end

end