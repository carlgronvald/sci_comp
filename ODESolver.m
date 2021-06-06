function [X,T] = ODESolver(x0, f, jac, h0, t0, t1, abstol, reltol, solver, parameters)
%ODESOLVER - Solves an ODE IVP using the given solver
% x0 - initial point
% f - handle for the ODE function, dx = f(t,x,parameters)
% jac - ODE function jacobian, only used for implicit methods.
% h0 - Starting step size. For fixed step size methods, this is the step
% size across the board.
% t0 - starting time.
% t1 - time we're approximating to
% abstol - Absolute tolerance of adaptive step size methods, not used by
% fixed step size methods
% reltol - Relative tolerance of adaptive step size methods, not used by
% fixed step size methods.
% solver - the solver you wish to use, see below.
% parameters - parameters for the ODE function and its jacobian
%
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

