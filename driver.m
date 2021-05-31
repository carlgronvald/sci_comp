%W = Wiener(0, 1, 1000, 2, 1, 0);
%Wsmall = zeros(2,100);
%for i=1:100
%    Wsmall(:,i) = W(:,i*10);
%end
%tic
%for i=0:100000
%    [X, T] = ExplicitEulerFixedStepSize([1.], @JacobianTest, 100, 0, 10, 0);
%end
%toc
%[X, T] = ExplicitEulerStepDoubling([1. 1.], @EquationTest, 0.1, 0, 1, [0.01 0.01], [0.01 0.01], 0)
%[X,T] = ImplicitEulerFixedStepSize([1.], @JacobianTest, 20, 0, 1, 1.0e-8, 100, 0)
%[X,T] = ImplicitEulerStepDoubling([1.], @JacobianTest, 0.01, 0, 1, 1.0e-8, 100, [0.001], [0.001], 0)
%[X, T] = RK4FixedStepSize([1.], @EquationTest, 20, 0, 1, 0)
%[X,T] = Dopri54AdaptiveStepSize([1.;2.0], @EquationTest, 0.1, 0, 2, [0.00001], [0.00001], 0);
%% Global and Local error, Explicit & Implicit Euler

max_steps = 1000;
local_errors = zeros(2,1);
global_errors = zeros(2,1);
hs = [];

for i=100:max_steps
    h = 5.0/i;
    hs = [hs;h];
    [X1,~] = ExplicitEulerFixedStepSize([1.0], @testf, i, 0, 5, 1);
    [X2,T] = ImplicitEulerFixedStepSize([1.0], @testf, @testjac, i, 0, 5, 1);
    d1 = diff(X1);
    d2 = diff(X2);
    dreal = diff(exp(T));
    local_exp = max(abs(d1-dreal));
    local_imp = max(abs(d2-dreal));
    local_errors = [local_errors [local_exp;local_imp]];
    global_errors = [global_errors [abs(X1(end)-exp(5));abs(X2(end)-exp(5))]];
    
    %TODO: SOMETHING VERY WEIRD IS HAPPENING WITH IMPLICIT
    %If we go small enough, it is all good.
end

plot(hs, global_errors(1:2,2:end));
legend(["Explicit", "Implicit"])
title("Global error of explicit and implicit euler at t=5, test equation")
xlabel("h")
ylabel("Global error")
shg


%plot(hs, local_errors(1:2,2:end));
%legend(["Explicit", "Implicit"])
%title("highest local error of explicit and implicit euler up to t=5, test equation")
%xlabel("h")
%ylabel("Maximum local error")
%shg


%[X, T] = ExplicitEulerStepDoubling([1. 1.], @testf, 0.1, 0, 1, [0.01 0.01], [0.01 0.01], 0)
%[X,T] = ImplicitEulerFixedStepSize([1.], @testf, @testjac, 20, 0, 1, 1.0e-8, 100, 0)
%[X,T] = ImplicitEulerStepDoubling([1.], @testf, @testjac, 0.01, 0, 1, 1.0e-8, 100, [0.001], [0.001], 0)
%[X, T] = RK4FixedStepSize([1.], @testf, 20, 0, 1, 0)
%[X, T] = RK4StepDoubling([1.], @testf, 0, 1, 0.1, 0.0001, 0.0001, 0);
%[X,T] = Dopri54AdaptiveStepSize([1.], @testf, 0.1, 0, 2, [0.00001], [0.00001], 0)
%[X,T] = ESDIRK23([1.], @testf, @testjac, 0, 1, 0.1, 0.001, 0.001, [])


%W = Wiener(0, 1, 1000, 2, 1, 0)
%Wsmall = zeros(2,100);
%for i=1:100
%    Wsmall(:,i) = W(:,i*10);
%end

%THIS SHOULD BE DONE IN A SMARTER WAY
%[X1, T] = SDEExplicitExplicitFixedStepSize(@testf, @testg, 1, 0, 1000, [1,2], W, 0);
%[X2, T] = SDEExplicitExplicitFixedStepSize(@testf, @testg, 1, 0, 99, [1,2], Wsmall, 0);
%Err = (X1(1:10:1000,:)-X2)