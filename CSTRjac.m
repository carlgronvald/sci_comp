function [jac] = CSTRjac(~, x, p)
% Jacobian of the 3D adiabatic continuous stirred tank reactor problem from the
%supplied paper
%x(1) = CA
%x(2) = CB
%x(3) = T
% Takes parameters F (flow rate), V (volume of tank), k0, EAonR (activation
% energy divided by the gas constant), CAin (input concentration of A, CBin
% (input concentration of B), and Tin (temperature of incoming flow)
F = p('F');
V = p('V');
beta = p('beta');
k0 = p('k0');
EaonR = p('EaonR');
%CAin = p('CAin');
%CBin = p('CBin');
%Tin = p('Tin');
T = x(3);
CA = x(1);
CB = x(2);

jac = [-F/V - CB*k0*exp(-EaonR/T), -CA*k0*exp(-EaonR/T), -CA*CB*k0*EaonR*exp(-EaonR/T)/T^2; ...
    -2*CB*k0*exp(-EaonR/T), -F/V - 2*CA*k0*exp(-EaonR/T), -2*CA*CB*k0*EaonR*exp(-EaonR/T)/T^2 ; ...
    beta*CB*k0*exp(-EaonR/T), beta*CA*k0*exp(-EaonR/T), -F/V + beta*CA*CB*k0*EaonR*exp(-EaonR/T)/T^2];
%jac = [ -fv-x(2)*k, -x(1)*k , EaonR*k*x(1)*x(2); ...
%        -2*x(2)*k, -fv-2*x(1)*k , EaonR*k*2*x(1)*x(2); ...
%        parameters('beta')*x(2)*k, parameters('beta')*x(1)*k, -fv-parameters('beta')*EaonR*k*x(1)*x(2)];
end

