function [dx] = CSTR1Df(t, x, p)
%CSTR1DF Summary of this function goes here
%   Detailed explanation goes here
F = p('F');
V = p('V');
beta = p('beta');
k0 = p('k0');
EaonR = p('EaonR');
CAin = p('CAin');
CBin = p('CBin');
Tin = p('Tin');
T = x;

dx = F*(Tin - T)/V + beta*k0*exp(-EaonR/T)*(CAin + (Tin - T)/beta)*(CBin + 2*(Tin - T)/beta);
end

