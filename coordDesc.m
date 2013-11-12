function [B, slackVars] = coordDesc(X, Y, lambda, gamma1, gamma2, coordFuncs, objFunc, Psi, Theta)

%coordFuncs is a cell array: the first entry contains the updates for B,
%and the additional entries contain any slack variable updates


tolerance = 1e-6;%Unsure

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
prevObj = 10e12;
numIT = 0;
while abs(currObj - prevObj) > tolerance
    numIT = numIT + 1;
    for i = 1:(numSlack+1)
        if i == 1
            B = feval(coordFuncs{i}, X, Y, B, Psi, Theta, slackVars, lambda, gamma1, gamma2);
            if any(isnan(B))
                keyboard;
            end
        else
            slackVars{i-1} = feval(coordFuncs{i},  X, Y, B, Psi, Theta, slackVars, lambda, gamma1, gamma2);
            if any(isnan(slackVars{i-1}))
                keyboard;
            end
        end
    end
    if numIT > 500
%         keyboard;
        disp(abs(currObj-prevObj))
    end
    prevObj = currObj;
    currObj = feval(objFunc, X, Y, B, slackVars, lambda, gamma1, gamma2, Psi, Theta);
    if currObj < 0
        keyboard;
    end
%     disp('meow')
% fprintf('Diff = %d\n', abs(currObj-prevObj));
% fprintf('Obj = %d\n', currObj);
end
% fprintf('Diff = %d\n', abs(currObj-prevObj));
fprintf('%d Iterations\n', numIT);
end