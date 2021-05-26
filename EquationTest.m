function [dxdt] = EquationTest(~, x, ~)
% EquationTest
% This is the differential equation dx/dt = x. Returns dx=x
% Doesn't use t and params for anything(Therefore, they're blank)
    dxdt = x;
end

