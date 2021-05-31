function [X, T] = RK4StepDoubling(x0, f, t0, t1, h0, abstol, reltol, params)
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
Butcher = [0 1/2 1/2 1 0; 0 1/2 0 0 1/6; 0 0 1/2 0 1/3; 0 0 0 1 1/3; 0 0 0 0 1/6]';

function x = KuttaStep(curx, curt, h)
    KuttaConstants(:,1) = curx;
    KuttaTimes = curt + h*Butcher(1:4,1);
    for s = 2:4
        KuttaConstants(:,s) = curx + h*Butcher(s, 2:s) * f(KuttaTimes(1:s-1), KuttaConstants(1, 1:s-1), params)';
    end
    x = curx + h*Butcher(5, 2:5) * f(KuttaTimes, KuttaConstants, params)';
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
            
            T = [T;t];
            X = [X;x];
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
