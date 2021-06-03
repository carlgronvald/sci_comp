function [dxdW] = CSTR1Dg(t,x,p)
%CSTR1DG Summary of this function goes here
%   Detailed explanation goes here
dxdW = p('F')/p('V') * p('sigma');
end

