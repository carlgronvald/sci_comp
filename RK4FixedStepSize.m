function [X, T] = RK4FixedStepSize(x0, f, h, t0, t1, params)
%RK4FIXEDSTEPSIZE - Fixed step size variant of the classical Runge Kutta
%method, a fourth order explicit RK method.
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

KuttaConstants = zeros(length(x0),4); %KuttaConstants are the internal Xs
% for each stage
KuttaTimes = zeros(1,4);
%Premultiply the butcher table with h, since it's fixed
hButchercs = h*[0 1/2 1/2 1];
hButcherbs = h*[1/6 1/3 1/3 1/6];
%Every stage only has a single term, so this is the terms in the butcher
%table that each stage uses.
hButcheras = h*[0 1/2 1/2 1];

for i = 1:steps
    if T(i)+h > t1
        h = t1-T(i);
    end
    KuttaConstants(:,1) = X(:,i);
    KuttaTimes = T(i) + hButchercs(1:4);
    for s = 2:4 %Calculate X for each stage
        KuttaConstants(:,s) = X(:,i) + hButcheras(s) * ...
            f(KuttaTimes(s-1), KuttaConstants(:, s-1), params);
    end
    
    %Sum the stages
    X(:,i+1) = X(:,i) +  hButcherbs(1) * f(KuttaTimes(1), KuttaConstants(:,1), params) ...
        + hButcherbs(2) * f(KuttaTimes(2), KuttaConstants(:,2), params) ...
        + hButcherbs(3) * f(KuttaTimes(3), KuttaConstants(:,3), params) ...
        + hButcherbs(4) * f(KuttaTimes(4), KuttaConstants(:,4), params);
    T(i+1) = T(i)+h;
end
X = X';
T = T';
