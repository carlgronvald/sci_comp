% 1 == Approximate test equation to t=5
% 2 == Oscillating euler CSTR approximation
mode = 1;

%% TEST EQUATION t=5
if mode == 1
%Here, we approximate the test equation up to t=5
solvers = ["ExplicitEulerFixedStepSize";"ExplicitEulerStepDoubling";"ImplicitEulerFixedStepSize"; ...
    "ImplicitEulerStepDoubling";"RK4FixedStepSize";"RK4StepDoubling";"DOPRI54";"ESDIRK23"];
%Below are the settings for the different solvers in order to achieve an
%error at t=5 of approximately 2.5.
%     EE     EEA     IE     IEA   RK4   RK4A  DO    ES
h0 = [0.0073,  0.01, 0.0072, 0.01, 1.007,    1,  1,  0.01];
tol =[0.01, 0.000065, 0.01,  0.000065, 0.01, 0.3, 0.5, 0.000062];
x0 = [1.0;1.0];
counts = zeros(8,1);
Jcounts = zeros(8,1);
times = zeros(8,1);
err = zeros(8,1);
steps = zeros(8,1);
%Count number of f and jacobian calls to compare efficiency.
global counter;
global Jcounter;
figure
hold on
for i=1:8
    counter = 0;
    Jcounter = 0;
    tic
    [X,T] = ODESolver(1.0,@testcounter, @testJcounter, h0(i), 0, 5, tol(i), tol(i), convertStringsToChars(solvers(i)), 1);
    times(i) = toc; %Wall time is used since the time is too small for CPU time to give reasonable results
    plot(T, X);
    counts(i) = counter;
    Jcounts(i) = Jcounter;
    err(i) = abs(X(end)-exp(5));
    steps(i) = length(T);
end
legend(solvers)
title("The test equation dx=x")
xlabel("t")
ylabel("x(t)")
disp([counts Jcounts times err])

end

%% CSTR1D explicit euler oscillating
if mode==2
% Construct a situation where Explicit Euler oscillates solving the 1D CSTR 
% problem, instead of finding the correct temperature
parameters = CSTRparameters();
x0 = 273.15;
parmcstr = @(t,x) CSTR1Df(t,x,parameters);
[X1,T1] = ExplicitEulerFixedStepSize(x0, @CSTR1Df, 50, 0, 400, parameters);
[Tcorrect, Xcorrect] = ode15s(parmcstr, [0 400], x0);

hold off
plot(T1, X1)
hold on
plot(Tcorrect, Xcorrect)
xlabel("t [s]")
ylabel("T [K]")
title("1D CSTR, temperature, explicit Euler, oscillations")
legend("fixed step", "ode15s")
end

%Counter functions to compare efficiency
function dx = testcounter(t,x,p)
    global counter;
    dx = testf(t,x,p);
    counter = counter+1;
end
function J = testJcounter(t,x,p)
     global Jcounter;
     J = testjac(t,x,p);
     Jcounter = Jcounter+1;
end
