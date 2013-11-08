function slackVar = fusedFusedSlack2Update(~, ~, B, Psi, ~, slackVars, lambda, gamma1, gamma2) %#ok<*INUSD>

[p, q] = size(B);

slackVar = zeros(q, p, p);
normalize = 0;
for k = 1:q
    for l = 1:p
        for m = 1:p
            slackVar(k, l, m) = Psi(l,m)*abs(B(l,k) - sign(Psi(l,m))*B(m, k));
            normalize = normalize + Psi(l,m)*abs(B(l,k) - sign(Psi(l,m))*B(m, k));
        end
    end
end

slackVar = slackVar/normalize;

end