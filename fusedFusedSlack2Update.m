function slackVar = fusedFusedSlack2Update(~, ~, B, Psi, ~, slackVars, lambda, gamma1, gamma2) %#ok<*INUSD>

[p, q] = size(B);

slackVar = slackVars{2};
normalize = 0;
for k = 1:q
    for l = 1:p
        for m = 1:p
            if (m~=l)
                slackVar(k, l, m) = Psi(l,m)*abs(B(l,k) - sign(Psi(l,m))*B(m, k));
                normalize = normalize + Psi(l,m)*abs(B(l,k) - sign(Psi(l,m))*B(m, k));
            else
                slackVar(k,l,m) = 1e-8;
                normalize = normalize + 1e-8;
            end
            if abs(slackVar(k,l,m)) <= 10^(-150)
%                 keyboard;
            end
        end
    end
end

slackVar = slackVar/normalize;

end