function [X,T] = ODESolver(x0, f, jac, h0, t0, t1, abstol, reltol, solver, parameters)
%ODESOLVER Summary of this function goes here
%   Detailed explanation goes here
% solvers: 'ExplicitEulerFixedStepSize', 'ExplicitEulerStepDoubling',
% 'ImplicitEulerFixedStepSize', 'ImplicitEulerStepDoubling',
% 'Dopri54AdaptiveStepSize', 'RK4FixedStepSize', 'RK4StepDoubling',
% 'ESDIRK23'

if strcmp('ExplicitEulerFixedStepSize', solver)
    %if abstol ~= 0 || reltol ~= 0
    %    disp("Explicit Euler Fixed Step Size does not need tolerances!")
    %end
    %if jac ~= 0
    %    disp("Explicit Euler Fixed Step Size does not need a jacobian!")
    %end
    [X,T] = ExplicitEulerFixedStepSize(x0, f, h0, t0, t1, parameters);
elseif strcmp('ExplicitEulerStepDoubling', solver)
    %if jac ~= 0
    %    disp("Explicit Euler Step Doubling does not need a jacobian!")
    %end
    [X,T] = ExplicitEulerStepDoubling(x0, f, h0, t0, t1, abstol, reltol, parameters);
elseif strcmp('ImplicitEulerFixedStepSize', solver)
    %if abstol ~= 0 || reltol ~= 0
    %    disp("Implicit Euler Fixed Step Size does not need tolerances!")
    %end
    [X,T] = ImplicitEulerFixedStepSize(x0, f, jac, h0, t0, t1, parameters);
elseif strcmp('ImplicitEulerStepDoubling', solver)
    [X,T] = ImplicitEulerStepDoubling(x0, f, jac, h0, t0, t1, abstol, reltol, parameters);
elseif strcmp('RK4FixedStepSize', solver)
    %if abstol ~= 0 || reltol ~= 0
    %    disp("RK4 Fixed Step Size does not need tolerances!")
    %end
    %if jac ~= 0
    %    disp("RK4 Fixed Step Size does not need a jacobian!")
    %end
    [X,T] = RK4FixedStepSize(x0, f, h0, t0, t1, parameters);
elseif strcmp('RK4StepDoubling', solver)
    %if jac ~= 0
    %    disp("RK4 Step Doubling does not need a jacobian!")
    %end
    [X,T] = RK4StepDoubling(x0, f, h0, t0, t1, abstol, reltol, parameters);
elseif strcmp('DOPRI54', solver)
    %if jac ~= 0
    %    disp("DOPRI54 does not need a jacobian!")
    %end
    [X,T] = Dopri54(x0, f, h0, t0, t1, abstol, reltol, parameters);
elseif strcmp('ESDIRK23', solver)
    [X,T] = ESDIRK23(x0, f, jac, h0, t0, t1, abstol, reltol, parameters);
else
    error(strcat("Unknown solver ", solver))
end




end

