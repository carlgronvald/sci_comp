function [dxdw] = vanderpolg(t,x,parameters)
%VANDERPOLG g function for the Van der Pol SDE with spread sigma 
dxdw = zeros(2,1);
dxdw(2) = parameters('sigma');
end

