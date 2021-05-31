W = Wiener(0, 1, 1000, 2, 1, 0);
Wsmall = zeros(2,100);
for i=1:100
    Wsmall(:,i) = W(:,i*10);
end
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

%THIS SHOULD BE DONE IN A SMARTER WAY
[X1, T] = SDEExplicitExplicitFixedStepSize(@testf, @testg, 1, 0, 1000, [1,2], W, 0);
[X2, T] = SDEExplicitExplicitFixedStepSize(@testf, @testg, 1, 0, 99, [1,2], Wsmall, 0);
Err = (X1(1:10:1000,:)-X2)