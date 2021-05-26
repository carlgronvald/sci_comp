function [X,T] = ExplicitEulerStepDoubling(x0, f, h0, t0, t1, abstol, reltol, params)
%EXPLICITEULERSTEPDOUBLING Computes explicit euler, adapting step size by
%using step doubling.
%   abstol = how large can the absolute error be (for each variable x)
%   reltol = how large can the relative error be (for each variable x)


epstol = 0.8; %epstol = what part of the maximal step we'll take next time TODO: BETTER DESCRIPTION
facmin = 0.1; %facmin = the smallest factor we'll allow multiplying h with in each step
facmax = 5.0; %facmax = the largest factor we'll allow multiplying h with in each step

h = h0;
t = t0;
x = x0;
X = x0;
T = t0;

while t < t1
    if t+h > t1
        h = t1-t;
    end
    
    AcceptStep = false;
    
    xdot = f(t, x, params);
    
    while ~AcceptStep %Keep trying until we find a step with sufficiently small error.
        x1 = x + h*xdot;
        
        hhalf = 0.5*h;
        thalfstep = t + hhalf;
        xhalfstep = x + xdot*hhalf;
        xdothalfstep = f(thalfstep, xhalfstep, params);
        x1doublestep = xhalfstep + xdothalfstep * hhalf;
        
        % Actual error
        e = x1doublestep-x1;
        % For each x variable, find the relation between either the
        % absolute tolerance or the relative tolerance times the estimated
        % correct value (whichever is best?? WHY BEST TODO). Use the worst of those as r
        % (kind of a how-bad-are-we but normalized estimate)
        r = max(abs(e)./max(abstol, x1doublestep.*reltol));
        
        AcceptStep = (r <= 1.0);
        if AcceptStep
            t = t+h;
            x = x1doublestep; %use the better estimate of the true value
            
            T = [T;t];
            X = [X;x];
        end
        
        % Calculate sqrt(epstol/r) - this is the 'largest step' (because of
        % 2nd order TODO). Then clamp between facmin and facmax, since we
        % don't change h by anymore than those factors, and multiply h
        % by them.
        h = max(facmin, min( sqrt(epstol/r), facmax))*h;
    end
end
X = X';
T = T';

