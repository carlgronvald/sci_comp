function [X, T] = ImplicitEulerStepDoubling(x0, f,jac, h0, t0, t1, abstol, reltol, params)

if size(x0,2) > 1
    error("x0 should be pased as column vector!")
end

epstol = 0.8; %epstol = what part of the maximal step we'll take next time TODO: BETTER DESCRIPTION
facmin = 0.1; %facmin = the smallest factor we'll allow multiplying h with in each step
facmax = 5.0; %facmax = the largest factor we'll allow multiplying h with in each step

newtonTolerance = 1.0e-8;
newtonMaxiterations = 100;

variable_count = size(x0);
variable_count = variable_count(2);

x = x0;
t = t0;
h = h0;
X = x;
T = t;
while t < t1
    if t+h > t1
        h = t1-t;
    end
    
    xdot = f(t, x, params);
    AcceptStep = false;
    
    
    while ~AcceptStep %Keep trying until we find a step with sufficiently small error.
        xguess = x+xdot*h;
        [x1, ~] =  NewtonsMethod(f, jac,  t, x, h, xguess, newtonTolerance, newtonMaxiterations, params);
        
        hhalf = 0.5*h;
        thalfstep = t + hhalf;
        [xhalfstep, ~] = NewtonsMethod(f, jac, t, x, hhalf, x + xdot*hhalf, newtonTolerance, newtonMaxiterations, params);
        xdothalfstep = f(thalfstep, xhalfstep, params);
        [x1doublestep, ~] = NewtonsMethod(f, jac, t, xhalfstep, hhalf, xhalfstep + xdothalfstep*hhalf, newtonTolerance, newtonMaxiterations, params);
        %TODO: CAN I USE THE dxdt FROM NEWTON'S METHOD??
        
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
            
            T = [T,t];
            X = [X,x];
        end
        
        % Calculate sqrt(epstol/r) - this is the 'largest step' (because of
        % 2nd order TODO). Then clamp between facmin and facmax, since we
        % don't change h by anymore than those factors, and multiply h
        % by them.
        h = max(facmin, min( sqrt(epstol/r), facmax))*h;
    end
end
while t < t1
end

T = T';
X = X';
