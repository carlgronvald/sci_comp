function [X, T] = SDEExplicitExplicitFixedStepSize(x0, f, g, h, t0, t1, W, params)
%SDEEXPLICITEXPLICITFIXEDSTEPSIZE Summary of this function goes here
%   Detailed explanation goes here
if size(x0,2) > 1
    error("x0 should be pased as column vector!")
end

steps = ceil((t1-t0)/h);
dimensions = length(x0);
X = zeros(dimensions, steps);
T = zeros(1, steps+1);

X(:,1) = x0;
T(1) = t0;
for k=1:steps
    if T(k)+h > t1
        h = t1-T(k);
    end
    dW = W(:,k+1)-W(:,k);
    dxdt = f(T(k), X(:,k), params);
    dxdw = g(T(k), X(:,k), params);
    X(:,k+1) = X(:,k) + dxdt*h + dxdw.*dW;
    T(k+1) = T(k)+h;
end
X = X';
T = T';