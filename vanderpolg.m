function [dxdw] = vanderpolg(t,x,parameters)
%VANDERPOLG Summary of this function goes here
%   Detailed explanation goes here
dxdw = zeros(2,1);
dxdw(1) = parameters('sigma1');
dxdw(2) = parameters('sigma2');
end

