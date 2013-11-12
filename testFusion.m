%% Initialization
clear;
clc;

% Functions to pass to coordinate descent
objFunc = @fusedObjective;
coordFuncs = {@fusedFusedUpdate; @fusedFusedSlack1Update; @fusedFusedSlack2Update; @fusedFusedSlack3Update};

%% Sanity check: Test the fused coordinate descent and compare resulting B

% Parameters for Synthetic Data
n = 100; %number of samples
p = 10; %number of features
q = 30; %number of outputs

density = 0.5; %density of trueB

[X, Y, Psi, Theta, trueB] = synthData(n, p, q, density);

figure
imagesc(trueB)
title('truth');
tic
[B, slackVar, lambda, gamma1, gamma2] = coordDescReg(X, Y, coordFuncs, objFunc, Psi, Theta);
toc
fprintf('Learned B:\n');
disp(B);
fprintf('True B:\n');
disp(trueB);

%% Experiment 1: How time and accuracy vary with n
% Parameters for Synthetic Data
p = 10;
q = 10;

density = 0.5;

N = 1000;
times = zeros(1, N/10);
accs = zeros(1, N/10);
t = 1;
for n = 10:10:N
    
    [X, Y, Psi, Theta, trueB] = synthData(n, p, q, density);
    
    time = cputime;
    [B, slackVar, lambda, gamma1, gamma2] = coordDescReg(X, Y, coordFuncs, objFunc, Psi, Theta);
    times(t) = cputime - time;
    accs(t) = norm(B - trueB);
    t = t+1;
end

save fusedFusedVaryN times accs;
figure;
plot(10:10:N, accs, 'r')
xlabel('Number of Samples');
ylabel('||B - trueB||');
title('Accuracy versus Number of Samples');
figure;
plot(10:10:N, times, 'b')
xlabel('Number of Samples');
ylabel('Runtime (s)');
title('Runtime versus Number of Samples');
%% Experiment 2: How time and accuracy vary with p

% Parameters for Synthetic Data
n = 100;
q = 10;

density = 0.5;

P = 1000;
times = zeros(1, P/10);
accs = zeros(1, P/10);
t = 1;
for p = 10:10:P
    
    [X, Y, Psi, Theta, trueB] = synthData(n, p, q, density);
    
    time = cputime;
    [B, slackVar, lambda, gamma1, gamma2] = coordDescReg(X, Y, coordFuncs, objFunc, Psi, Theta);
    times(t) = cputime - time;
    accs(t) = norm(B - trueB);
    t = t+1;
end

save fusedFusedVaryP times accs;
figure;
plot(10:10:P, accs, 'r')
xlabel('Number of Features');
ylabel('||B - trueB||');
title('Accuracy versus Number of Features');
figure;
plot(10:10:P, times, 'b')
xlabel('Number of Features');
ylabel('Runtime (s)');
title('Runtime versus Number of Features');
%% Experiment 3: How time and accuracy vary with q

% Parameters for Synthetic Data
n = 100;
p = 10;

density = 0.5;

Q = 1000;
times = zeros(1, Q/10);
accs = zeros(1, Q/10);
t = 1;
for q = 10:10:Q
    
    [X, Y, Psi, Theta, trueB] = synthData(n, p, q, density);
    
    time = cputime;
    [B, slackVar, lambda, gamma1, gamma2] = coordDescReg(X, Y, coordFuncs, objFunc, Psi, Theta);
    times(t) = cputime - time;
    accs(t) = norm(B - trueB);
    t = t+1;
end

save fusedFusedVaryP times accs;
figure;
plot(10:10:P, accs, 'r')
xlabel('Number of Outputs');
ylabel('||B - trueB||');
title('Accuracy versus Number of Outputs');
figure;
plot(10:10:Q, times, 'b')
xlabel('Number of Outputs');
ylabel('Runtime (s)');
title('Runtime versus Number of Outputs');
%% Experiment 4: How time and accuracy vary with density

% Parameters for Synthetic Data
n = 100;
p = 10;
q = 10;


D = 1;
times = zeros(1, 100);
accs = zeros(1, 100);
t = 1;
for density = 0.1:0.1:D
    
    [X, Y, Psi, Theta, trueB] = synthData(n, p, q, density);
    
    time = cputime;
    [B, slackVar, lambda, gamma1, gamma2] = coordDescReg(X, Y, coordFuncs, objFunc, Psi, Theta);
    times(t) = cputime - time;
    accs(t) = norm(B - trueB);
    t = t+1;
end

save fusedFusedVaryD times accs;
figure;
plot(0.1:0.1:1, accs, 'r')
xlabel('Density');
ylabel('||B - trueB||');
title('Accuracy versus Density');
figure;
plot(0.1:0.1:1, times, 'b')
xlabel('Density');
ylabel('Runtime (s)');
title('Runtime versus Density');
