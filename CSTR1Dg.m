function [dxdW] = CSTR1Dg(t,x,p)
%CSTR1DG Returns the derivative of the temperature with regard to the
%derivative of the Wiener process. Takes a parameter sigma for the spread
%of the temperature.
dxdW = p('F')/p('V') * p('sigma');
end

