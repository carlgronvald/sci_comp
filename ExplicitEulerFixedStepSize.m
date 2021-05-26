function [X, T] = ExplicitEulerFixedStepSize(x0, f, steps, t0, t1, params)
h = (t1-t0)/steps;
variable_count = size(x0);
variable_count = variable_count(2);

X = zeros(variable_count, steps+1);
X(:,1) = x0;
T = zeros(1, steps+1);
T(1) = t0;


for i = 1:steps
    X(:,i+1) = X(:,i) + h*f(T(i), X(:,i), params);
    T(i+1) = T(i)+h;
end
X = X';
T = T';
