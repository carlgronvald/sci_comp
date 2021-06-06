function [dx] = testf(t, x, parameters)
%TESTF ODE function for the ODE dx= lambda x
dx = x*parameters;
end

