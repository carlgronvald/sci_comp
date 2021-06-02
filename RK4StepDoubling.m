function [X, T] = RK4StepDoubling(x0, f, h0, t0, t1, abstol, reltol, params)
if size(x0,2) > 1
    error("x0 should be pased as column vector!")
end
h = h0;
X = [x0];
T = t0;
t = t0;
x = x0;

epstol = 0.8; %epstol = what part of the maximal step we'll take next time TODO: BETTER DESCRIPTION
facmin = 0.1; %facmin = the smallest factor we'll allow multiplying h with in each step
facmax = 5.0; %facmax = the largest factor we'll allow multiplying h with in each step

KuttaConstants = zeros(length(x0),4);
KuttaTimes = zeros(1,4);

Butcheras = [0 1/2 1/2 1];
Butchercs = [0 1/2 1/2 1];
Butcherbs = [1/6 1/3 1/3 1/6];
function x = KuttaStep(curx, curt, h)
    
    KuttaConstants(:,1) = curx;
    KuttaTimes = curt + h*Butchercs;
    for s = 2:4
        KuttaConstants(:,s) = curx + h*Butcheras(s) * f(KuttaTimes(s-1), KuttaConstants(:, s-1), params);
    end
    
   	x = curx +  h*Butcherbs(1) * f(KuttaTimes(1), KuttaConstants(:,1), params) ...
        + h*Butcherbs(2) * f(KuttaTimes(2), KuttaConstants(:,2), params) ...
        + h*Butcherbs(3) * f(KuttaTimes(3), KuttaConstants(:,3), params) ...
        + h*Butcherbs(4) * f(KuttaTimes(4), KuttaConstants(:,4), params);
end

n=1;
while t<t1
    if t+h > t1
        h = t1-t;
    end
    AcceptStep = false;
    while ~AcceptStep
        large_step = KuttaStep(x, t, h);
        half_small_step = KuttaStep(x,t,h/2);
        small_step = KuttaStep(half_small_step,t+h/2,h/2);
        
        e = large_step-small_step;
        r = max(abs(e)./max(abstol, small_step.*reltol));
        
        AcceptStep = (r <= 1.0);
        if AcceptStep
            t = t+h;
            x = small_step; %use the better estimate of the true value
            
            T = [T,t];
            X = [X,x];
            n = n+1;
        end
        
        % Use asymptotic error control, our integrator is of order 4, so we
        % get ^(1/(p+1)) = ^(1/5)
        h = max(facmin, min( (epstol/r)^(1/5), facmax))*h;
    end
end
X = X';
T = T';

end
