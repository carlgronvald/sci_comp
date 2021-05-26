function [X,T] = SDEImplicitExplicitFixedStepSize(fJac, g, t1, t0, steps, x0, W, params)
%SDEEXPLICITIMPLICITFIXEDSTEPSIZE Summary of this function goes here
%   Detailed explanation goes here
%TODO: TOLERANCE AND MAXIT AS PARAMETERS
tol = 1.0e-8;
maxit = 100;

h = (t1-t0)/steps;
dimensions = length(x0);
X = zeros(dimensions, steps);
T = zeros(1, steps+1);


X(:,1) = x0;
T(1) = t0;
[dxdt,~] = fJac(T(0), X(:,0), params);
for k=1:steps
    dxdW = g(T(k), X(:,k), params);
    dW = W(:,k+1)-W(:,k);
    psi = X(:,k) + dxdW*dW;
    xguess = psi + dxdt*h;
    [X(:, k+1), dxdt] = NewtonsMethod(fJac, T(:,k+1), h, psi, xguess, tol, maxit, params);
end



