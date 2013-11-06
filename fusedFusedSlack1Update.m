function slackVar = fusedFusedSlack1Update(~, ~, B, ~, ~, slackVars, lambda, gamma1, gamma2) %#ok<*INUSD>
slackVar = abs(B)./sum(sum(abs(B)));

end