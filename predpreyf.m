function [dx] = predpreyf(t, x, params)
%PREDPREYF ODE function for the predator-prey problem.
dx = zeros(2,1);
dx(1) = params('a')*(1-x(2))*x(1);
dx(2) = -params('b')*(1-x(1))*x(2);
end

