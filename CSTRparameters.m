function [parameters] = CSTRparameters
%CSTRPARAMETERS Summary of this function goes here
%   Detailed explanation goes here

k = {'F', 'V',   'k0',      'EaonR', 'beta',           'CAin', 'CBin', 'Tin'};
v = {0.400/60, 0.105, exp(24.6), 8500,     560/(1.0*4.186), 0.8,    1.2,    273.65}; %TODO: UNITS

parameters = containers.Map(k, v);

end

