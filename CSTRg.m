function [dxdW] = CSTRg(~, ~, parameters)
%CSTRG Stochastic part of CSTR problem.
%sigma is the standard deviation of the temperature.

dxdW = [0;0; parameters('F')/parameters('V') * parameters('sigma')];
end

