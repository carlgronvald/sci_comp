function [parameters] = CSTRparameters(ek,ev)
%CSTRPARAMETERS Creates a default set of parameters for the CSTR problem

k = {'F', 'V',   'k0',      'EaonR', 'beta',           'CAin', 'CBin', 'Tin'};
v = {0.400/60, 0.105, exp(24.6), 8500,     560/(1.0*4.186), 0.8,    1.2,    273.65}; %TODO: UNITS
if(nargin == 2)
    k = [k,ek];
    v = [v,ev];
end

parameters = containers.Map(k, v);

end

