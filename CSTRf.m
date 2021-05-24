function [dx, jac] = CSTRf(x, parameters)
%CSTR The adiabatic continuous stirred tank reactor problem from the
%supplied paper
%x(1) = CA
%x(2) = CB
%x(3) = T
% Takes parameters F (flow rate), V (volume of tank), k0, EAonR (activation
% energy divided by the gas constant), CAin (input concentration of A, CBin
% (input concentration of B), and Tin (temperature of incoming flow)

fv = parameters('F')/parameters('V');
k0 = parameters('k0');
EaonR = parameters('Ea');

k = k0*exp(-EaonR/x(3));
r = x(1)*x(2)*k;
dx = [ fv*(parameters('CAin') - x(1)); ...
        fv*(parameters('CBin')-x(2)); ...
        fv*(parameters('Tin')-x(3))];
end

