function [X, T] = RK4FixedStepSize(x0, f, steps, t0, t1, params)
if size(x0,2) > 1
    error("x0 should be pased as column vector!")
end

h = (t1-t0)/steps;

X = zeros(length(x0), steps+1);
X(:,1) = x0;
T = zeros(1, steps+1);
T(1) = t0;

KuttaConstants = zeros(variable_count,4);
KuttaTimes = zeros(1,4);
hButcher = h*[0 1/2 1/2 1 0; 0 1/2 0 0 1/6; 0 0 1/2 0 1/3; 0 0 0 1 1/3; 0 0 0 0 1/6]';

for i = 1:steps
    KuttaConstants(:,1) = X(:,i);
    KuttaTimes = T(i) + hButcher(1:4, 1);
    for s = 2:4
        KuttaConstants(:,s) = X(:,i) + hButcher(s, 2:s) * f(KuttaTimes(1:s-1),KuttaConstants(1,1:s-1),params)';
    end
    
    X(:,i+1) = X(:,i) +  hButcher(5,2:5) * f(KuttaTimes, KuttaConstants, params)';
    T(i+1) = T(i)+h;
end
X = X';
T = T';
