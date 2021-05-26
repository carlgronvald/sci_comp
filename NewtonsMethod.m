function [x, xdot] = NewtonsMethod(fjac, told, Rterm, dt, xguess, tolerance, maxiterations, params)
%NEWTONSMETHODODE Does up to maxiterations rounds of newtons method to find
% the next step of the implicit within given tolerance.
i = 0;
t = told + dt;
x = xguess;
[xdot, J] = fjac(t, x, params);
R = x - xdot*dt - Rterm;
I = eye(length(Rterm));
while i < maxiterations && max(abs(R)) > tolerance %Iteratively improve guess using newton's method.
    i = i+ 1;
    dRdx = I - J * dt;
    dx = dRdx\R;
    x = x - dx;
    [xdot, J] = fjac(t, x, params);
    R = x - dt*xdot - Rterm;
end

