function [B, slackVars] = coordDesc(X, Y, lambda, gamma1, gamma2, coordFuncs, objFunc, Psi, Theta)

%coordFuncs is a cell array: the first entry contains the updates for B,
%and the additional entries contain any slack variable updates

tolerance = 1e-4;%Unsure

p = size(X,2);
q = size(Y,2);
numSlack = length(coordFuncs)-1;

B = rand(p, q);
slackVars = cell(1, numSlack);
%This is specific to fused-fused penalty
slackVars{1} = rand(p,q);
slackVars{2} = rand(q,p,p);
slackVars{3} = rand(p,q,q);

currObj = feval(objFunc, X, Y, B, slackVars, lambda, gamma1, gamma2, Psi, Theta);
prevObj = 10e10;
while abs(currObj - prevObj) > tolerance
    for i = 1:(numSlack+1)
        if i == 1
            B = feval(coordFuncs{i}, X, Y, B, Psi, Theta, slackVars, lambda, gamma1, gamma2);
        else
            slackVars{i-1} = feval(coordFuncs{i},  X, Y, B, Psi, Theta, slackVars, lambda, gamma1, gamma2);
        end
    end
    
    prevObj = currObj;
    currObj = feval(objFunc, X, Y, B, slackVars, lambda, gamma1, gamma2, Psi, Theta);
end
end