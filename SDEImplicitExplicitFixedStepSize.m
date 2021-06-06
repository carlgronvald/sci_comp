function [X,T] = SDEImplicitExplicitFixedStepSize(x0, f, jac, g, h, t0, t1, W, params)
%SDEEXPLICITIMPLICITFIXEDSTEPSIZE Solves an SDE using the Implicit-Explicit
%method, meaning it is explicit in W and implicit in t

if size(x0,2) > 1
    error("x0 should be pased as column vector!")
end

steps = ceil((t1-t0)/h);
%Instead of using the exact passed step size, we create a step size that
%takes the same number of steps but evenly, instead of cutting off the last
%step
h = (t1-t0)/steps;
tol = 1.0e-8;
maxit = 100;

dimensions = length(x0);
X = zeros(dimensions, steps);
T = zeros(1, steps+1);


X(:,1) = x0;
T(1) = t0;
dxdt = f(T(1), X(:,1), params);
for k=1:steps
    %calculate first the explicit contribution from dxdw
    dxdW = g(T(k), X(:,k), params);
    dW = W(:,k+1)-W(:,k);
    %add that to the residual
    psi = X(:,k) + dxdW.*dW;
    %Then solve the implicit equation x_n+1 = x_n + dw * g(x_n) + dt*f(x_n+1)
    xguess = psi + dxdt*h;
    T(k+1) = T(k)+h;
    [X(:, k+1), dxdt] = NewtonsMethod(f, jac, T(k), psi, h, xguess, tol, maxit, params);
end

T = T';
X = X';

