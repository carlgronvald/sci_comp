function [J] = vanderpoljac(t, x, parameters)
%VANDERPOLF Summary of this function goes here
%   Detailed explanation goes here
J = zeros(2,2);
J(1,1) = 0;
J(1,2) = 1;
J(2,1) = -2*parameters('mu')*x(2)*x(1) - 1;
J(2,2) = parameters('mu')*(1-x(1)^2);
end

