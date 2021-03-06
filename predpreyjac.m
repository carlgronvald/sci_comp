function [J] = predpreyjac(t, x, params)
%PREDPREYF Jacobian for the predator-prey problem.
J = zeros(2,2);
J(1,1) = params('a')*(1-x(2));
J(1,2) = -params('a')*x(1);
J(2,1) = params('b')*x(2);
J(2,2) = -params('b')*(1-x(1));
end

