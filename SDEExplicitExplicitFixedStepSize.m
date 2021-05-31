function [X, T] = SDEExplicitExplicitFixedStepSize(f, g, t1, t0, steps, x0, W, params)
%SDEEXPLICITEXPLICITFIXEDSTEPSIZE Summary of this function goes here
%   Detailed explanation goes here
if size(x0,2) > 1
    error("x0 should be pased as column vector!")
end

h = (t1-t0)/steps;
dimensions = length(x0);
X = zeros(dimensions, steps);
T = zeros(1, steps+1);

X(:,1) = x0;
T(1) = t0;
for k=1:steps
    dW = W(:,k+1)-W(:,k);
    dxdt = f(T(k), X(:,k), params);
    dxdw = g(T(k), X(:,k), params);
    X(:,k+1) = X(:,k) + dxdt*h + dxdw*dW;
    T(k+1) = T(k)+h;
end
X = X';
T = T';