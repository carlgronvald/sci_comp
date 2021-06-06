function [dx] = CSTR1Df(t, x, p)
%CSTR1DF Returns the derivative of the temperature in the one dimensional
%CSTR problem
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

