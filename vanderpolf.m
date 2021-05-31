function [dx] = vanderpolf(t, x, parameters)
%VANDERPOLF Summary of this function goes here
%   Detailed explanation goes here
dx = zeros(2, 1);
dx(1) = x(2);
dx(2) = parameters('mu') *( 1- x(1)^2) * x(2) - x(1);
end

