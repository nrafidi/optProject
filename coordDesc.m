function [B, slackVars, numIT] = coordDesc(X, Y, lambda, gamma1, gamma2, coordFuncs, objFunc, Psi, Theta, Bstart)

%coordFuncs is a cell array: the first entry contains the updates for B,
%and the additional entries contain any slack variable updates


tolerance = 1e-6;%Unsure

p = size(X,2);
q = size(Y,2);
numSlack = length(coordFuncs)-1;

B = Bstart;%rand(p, q);
slackVars = cell(1, numSlack);
%This is specific to fused-fused penalty
slackVars{1} = rand(p,q);
slackVars{2} = rand(q,p,p);
slackVars{3} = rand(p,q,q);

currObj = feval(objFunc, X, Y, B, slackVars, lambda, gamma1, gamma2, Psi, Theta);
prevObj = 10e20;
numIT = 0;
while abs(currObj - prevObj) > tolerance
    tic
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
    if numIT > 100
        %         keyboard;
        break;
    end
    prevObj = currObj;
    currObj = feval(objFunc, X, Y, B, slackVars, lambda, gamma1, gamma2, Psi, Theta);
    if currObj < 0
        keyboard;
    end
    if currObj - prevObj > 0 && numIT > 1
        %         keyboard;
        break;
    end
    toc
    %     disp('meow')
    % fprintf('Diff = %d\n', abs(currObj-prevObj));
    %     if gamma1 ~= 0
    %         fprintf('Obj = %d\n', currObj);
    %     end
end
% fprintf('Diff = %d\n', abs(currObj-prevObj));
% fprintf('%d Iterations\n', numIT);
end