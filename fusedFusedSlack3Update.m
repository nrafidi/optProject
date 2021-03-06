function slackVar = fusedFusedSlack3Update(~, ~, B, ~, Theta, slackVars, lambda, gamma1, gamma2) %#ok<*INUSD>

[p, q] = size(B);

slackVar = slackVars{3};
normalize = 0;
for j = 1:p
    for r = 1:q
        for s = 1:q
            if r ~=s
                slackVar(j, r, s) = abs(Theta(r,s))*abs(B(j,r) - sign(Theta(r,s))*B(j, s));
                normalize = normalize + abs(Theta(r,s))*abs(B(j,r) - sign(Theta(r,s))*B(j, s));
            else
                slackVar(j,r,s) = 1;%1e-8;
                normalize = normalize + 1;%1e-8;
            end
            if slackVar(j,r,s) <= 10^(-150)
                slackVar(j,r,s) = 10^(-100);
%                 keyboard;
            end
        end
    end
end

slackVar = slackVar/normalize;

end