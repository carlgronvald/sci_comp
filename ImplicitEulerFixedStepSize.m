function [X, T] = ImplicitEulerFixedStepSize(x0, f, jac, h, t0, t1, params)

if size(x0,2) > 1
    error("x0 should be pased as column vector!")
end
steps = ceil((t1-t0)/h);
newtonTolerance = 0.0000000000001;
maxiterations = 100;

newtonTolerance = 1.0e-8;
newtonMaxiterations = 100;


X = zeros(length(x0), steps+1);
X(:,1) = x0;
T = zeros(1, steps+1);
T(1) = t0;

for i = 1:steps
    if T(i)+h > t1
        h = t1-T(i);
    end
    xdot = f(T(i), X(:,i), params);
    T(i+1) = T(i)+h;
    xguess = X(:,i)+xdot*h;
    [X(:,i+1), ~] = NewtonsMethod(f, jac,  T(:,i), X(:,i), h, xguess, newtonTolerance, newtonMaxiterations, params);
end

T = T';
X = X';
