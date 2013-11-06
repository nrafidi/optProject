function [X, Y, Psi, Theta, trueB] = synthData(n, p, q, density)

trueB = sprandn(p, q, density);

Psi = iwishrnd(eye(p), p);
Theta = iwishrnd(eye(q), q);

cholPsi = chol(Psi, 'lower');
cholTheta = chol(Theta);

trueB = cholPsi*trueB*cholTheta;

X = mvnrnd(zeros(n, p), Psi);

Y = X*trueB;


end