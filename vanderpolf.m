function [dx] = vanderpolf(t, x, parameters)
%VANDERPOLF ODE function of the Van der Pol problem
dx = zeros(2, 1);
dx(1) = x(2);
dx(2) = parameters('mu') *( 1- x(1)^2) * x(2) - x(1);
end

