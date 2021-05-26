function [dxdt, jac] = JacobianLinear(~, x, ~)
%JACOBIANTEST Summary of this function goes here
%   Detailed explanation goes here
    dxdt = 5;
    shape = size(x);
    jac = zeros(shape(2), shape(2));
end

