function [X, T] = ExplicitEulerFixedStepSize(x0, f, h, t0, t1, params)

if size(x0,2) > 1
    error("x0 should be pased as column vector!")
end
steps = ceil((t1-t0)/h);
X = zeros(length(x0), steps+1);
X(:,1) = x0;
T = zeros(1, steps+1);
T(1) = t0;


for i = 1:steps
    if T(i)+h > t1
        h = t1-T(i);
    end
    X(:,i+1) = X(:,i) + h*f(T(i), X(:,i), params);
    T(i+1) = T(i)+h;
end
X = X';
T = T';
