function [dx] = predpreyf(t, x, params)
%PREDPREYF Summary of this function goes here
%   Detailed explanation goes here
dx = zeros(2,1);
dx(1) = params('a')*(1-x(2))*x(1);
dx(2) = -params('b')*(1-x(1))*x(2);
end

