function obj = fusedObjective(X, Y, B, slackVar, lambda, gamma1, gamma2, Psi, Theta)

obj = trace(Y - X*B) + lambda*sum(sum((B.^2)./slackVar{1})); %Not actually sure how to code the rest

end