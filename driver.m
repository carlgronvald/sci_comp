%tic
%for i=0:100000
%    [X, T] = ExplicitEulerFixedStepSize([1.], @JacobianTest, 100, 0, 10, 0);
%end
%toc
%[X, T] = ExplicitEulerStepDoubling([1. 1.], @EquationTest, 0.1, 0, 1, [0.01 0.01], [0.01 0.01], 0)
%[X,T] = ImplicitEulerFixedStepSize([1.], @JacobianTest, 20, 0, 1, 1.0e-8, 100, 0)
%[X,T] = ImplicitEulerStepDoubling([1.], @JacobianTest, 0.01, 0, 1, 1.0e-8, 100, [0.001], [0.001], 0)
%[X, T] = RK4FixedStepSize([1.], @EquationTest, 20, 0, 1, 0)
[X,T] = Dopri54AdaptiveStepSize([1.;2.0], @EquationTest, 0.1, 0, 2, [0.00001], [0.00001], 0);

%[X, T] = 