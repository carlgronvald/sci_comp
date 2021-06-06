function [X, T] = SDEExplicitExplicitFixedStepSize(x0, f, g, h, t0, t1, W, params)
%SDEEXPLICITEXPLICITFIXEDSTEPSIZE Solves an SDE using the Explicit-Explicit
%method, meaning it is explicit both in W and in t
if size(x0,2) > 1
    error("x0 should be pased as column vector!")
end

steps = ceil((t1-t0)/h);
%Instead of using the exact passed step size, we create a step size that
%takes the same number of steps but evenly, instead of cutting off the last
%step
h = (t1-t0)/steps;
dimensions = length(x0);
X = zeros(dimensions, steps);
T = zeros(1, steps+1);

X(:,1) = x0;
T(1) = t0;
for k=1:steps
    dW = W(:,k+1)-W(:,k);
    %calculate first the explicit contribution from dxdt, then from dxdw
    dxdt = f(T(k), X(:,k), params);
    dxdw = g(T(k), X(:,k), params);
    %Sum them
    X(:,k+1) = X(:,k) + dxdt*h + dxdw.*dW;
    T(k+1) = T(k)+h;
end
X = X';
T = T';