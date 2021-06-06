function [X, T] = ExplicitEulerFixedStepSize(x0, f, h, t0, t1, params)
%EXPLICITEULERFIXEDSTEPSIZE Solves a given ODE using the explicit euler
%method with fixed step size
if size(x0,2) > 1
    error("x0 should be pased as column vector!")
end
steps = ceil((t1-t0)/h);
%Instead of using the exact passed step size, we create a step size that
%takes the same number of steps but evenly, instead of cutting off the last
%step
h = (t1-t0)/steps;
X = zeros(length(x0), steps+1);
X(:,1) = x0;
T = zeros(1, steps+1);
T(1) = t0;


for i = 1:steps
    %Advance using the formula x_n+1 = x_n + hf
    X(:,i+1) = X(:,i) + h*f(T(i), X(:,i), params);
    T(i+1) = T(i)+h;
end
X = X';
T = T';
