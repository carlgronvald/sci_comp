function [jac] = testjac(t, x, parameters)
%TESTJAC Jacobian of the equation dx=lambda x
jac = eye(size(x))*parameters;

end

