subjects = 'ABCD';
load ../brainData/semantic_data.mat;

%% Coordinate Descent Code - Training
objFunc = @fusedObjective;
coordFuncs = {@fusedFusedUpdate; @fusedFusedSlack1Update; @fusedFusedSlack2Update; @fusedFusedSlack3Update};
paramOrder = [1 0 0 0 1]; %merrrh


betas = cell(length(subjects), 1);
for s = 1:length(subjects)
    load(['../brainData/' subjects(s) '_data.mat']);
    betas{s} = coordDescReg(X, Y, coordFuncs, objFunc, Psi, Theta, paramOrder);
    fprintf('Subject %s complete\n', subjects(s));
end

save ../brainData/betasCoord.mat betas;

%% Coordinate Descent Code - Testing and Plotting

if ~exist('betas', 'var')
    load ../brainData/beasCoord.mat;
end

N = size(Y,1);
S = length(subjects);

scores = zeros(S, S);

for s = 1:S
    subplot(S/2, S/2, s);
    imagesc(betas{s});
    colorbar;
    xlabel('Semantic Features');
    ylabel('MEG Sources');
    title(subjects(s));
    for ss = 1:S
        if ss ~= s
            load(['../brainData/' subjects(ss) '_data.mat']);
            Yhat = X*betas{s};
            for n = 1:N
                d = pdist([Yhat(n,:); Y], 'cosine');
                [~,ind] = sort(d(1:N));
                ranks(n) = find(n==ind);
            end
            scores(s, ss) = mean(1 - (ranks - 1)./(N -1));
        end
    end
end
suptitle('Weight maps for individual subjects');
figure;
imagesc(scores);
colorbar;
xlabel('Testing Subject');
ylabel('Training Subject');
title('Rank Accuracy Scores Across Subjects');