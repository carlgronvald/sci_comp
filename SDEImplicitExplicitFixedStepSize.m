function [X,T] = SDEImplicitExplicitFixedStepSize(x0, f, jac, g, h, t0, t1, W, params)
%SDEEXPLICITIMPLICITFIXEDSTEPSIZE Summary of this function goes here
%   Detailed explanation goes here
%TODO: TOLERANCE AND MAXIT AS PARAMETERS

if size(x0,2) > 1
    error("x0 should be pased as column vector!")
end

steps = ceil((t1-t0)/h);
tol = 1.0e-8;
maxit = 100;

dimensions = length(x0);
X = zeros(dimensions, steps);
T = zeros(1, steps+1);


X(:,1) = x0;
T(1) = t0;
dxdt = f(T(1), X(:,1), params);
for k=1:steps
    if T(k)+h > t1
        h = t1-T(k);
    end
    dxdW = g(T(k), X(:,k), params);
    dW = W(:,k+1)-W(:,k);
    psi = X(:,k) + dxdW.*dW;
    xguess = psi + dxdt*h;
    T(k+1) = T(k)+h;
    [X(:, k+1), dxdt] = NewtonsMethod(f, jac, T(k), psi, h, xguess, tol, maxit, params);
end

T = T';
X = X';

