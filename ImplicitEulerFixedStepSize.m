function [X, T] = ImplicitEulerFixedStepSize(x0, f, jac, steps, t0, t1, newtonTolerance, maxiterations, params)
h = (t1-t0)/steps;

variable_count = size(x0);
variable_count = variable_count(2);

X = zeros(variable_count, steps+1);
X(:,1) = x0;
T = zeros(1, steps+1);
T(1) = t0;


for i = 1:steps
    xdot = f(T(i), X(:,i), params);
    T(i+1) = T(i)+h;
    xguess = X(:,i)+xdot*h;
    [X(:,i+1), ~] = NewtonsMethod(f, jac,  T(:,i), X(:,i), h, xguess, newtonTolerance, maxiterations, params);
end

T = T';
X = X';
