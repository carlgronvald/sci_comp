% 1 == Van der Pol RK4
% 2 == CSTR 3D RK4
% 3 == CSTR 1D RK4
% 4 == Order RK4
% 5 == RK4 Stability
mode = 5;

%% Van der Pol
if mode == 1
%Use RK4 both with step doubling and fixed step size on the Van der Pol
%problem, and compare to ode15s
%First, mu=1.5
parameters = CreateParams('mu', 1.5);
x0 = [1.0;1.0];
vanmu1p5 = @(t,x) vanderpolf(t,x,parameters);
[X1,T1] = RK4FixedStepSize(x0, @vanderpolf, 0.5, 0, 40, parameters);
[X2,T2] = RK4StepDoubling(x0, @vanderpolf, 0.2, 0, 40, 0.01, 0.01, parameters);
[Tcorrect, Xcorrect] = ode15s(vanmu1p5, [0 40], [1.0;1.0]);
hold off
plot(T1, X1(:,1))
hold on
plot(T2, X2(:,1))
plot(Tcorrect, Xcorrect(:,1))
title("Van der Pol, mu=1.5, Classical Runge Kutta")
xlabel("t")
ylabel("x(1)")
legend("fixed step", "adaptive step", "ode15s")

%Then, mu=15
figure 
parameters = CreateParams('mu', 15);
x0 = [1.0;1.0];
vanmu1p5 = @(t,x) vanderpolf(t,x,parameters);
[X1,T1] = RK4FixedStepSize(x0, @vanderpolf, 0.05, 0, 40, parameters);
[X2,T2] = RK4StepDoubling(x0, @vanderpolf, 0.2, 0, 40, 0.01, 0.01, parameters);
[Tcorrect, Xcorrect] = ode15s(vanmu1p5, [0 40], [1.0;1.0]);
hold off
plot(T1, X1(:,1))
hold on
plot(T2, X2(:,1))
plot(Tcorrect, Xcorrect(:,1))
title("Van der Pol, mu=15, Classical Runge Kutta")
xlabel("t")
ylabel("x(1)")
legend("fixed step", "adaptive step", "ode15s")
end
%% Adiabatic CSTR 3D
if mode == 2
%Use RK4 with fixed and adaptive step size on the adiabatic CSTR 3D problem
parameters = CSTRparameters();
x0 = CSTRx0(parameters);
parmcstr = @(t,x) CSTRf(t,x,parameters);
[X1,T1] = RK4FixedStepSize(x0, @CSTRf, 10, 0, 200, parameters);
[X2,T2] = RK4StepDoubling(x0, @CSTRf, 1, 0, 200, 0.01, 0.01, parameters);
[Tcorrect, Xcorrect] = ode15s(parmcstr, [0 200], x0);

hold off
plot(T1, X1(:,3))
hold on
plot(T2, X2(:,3))
plot(Tcorrect, Xcorrect(:,3))
xlabel("t [s]")
ylabel("T [K]")
title("3D CSTR, temperature, Classical Runge Kutta")
xlabel("t [s]")
ylabel("T [K]")
legend("fixed step", "step doubling", "ode15s")

end
%% Adiabatic CSTR 1D
if mode == 3
%RK4 with and without adaptive step size on CSTR 1D problem
parameters = CSTRparameters();
x0 = 273.15;
parmcstr = @(t,x) CSTR1Df(t,x,parameters);
[X1,T1] = RK4FixedStepSize(x0, @CSTR1Df, 10, 0, 200, parameters);
[X2,T2] = RK4StepDoubling(x0, @CSTR1Df, 1, 0, 200, 0.01, 0.01, parameters);
[Tcorrect, Xcorrect] = ode15s(parmcstr, [0 200], x0);

hold off
plot(T1, X1)
hold on
plot(T2, X2)
plot(Tcorrect, Xcorrect)
xlabel("t [s]")
ylabel("T [K]")
title("1D CSTR, temperature, RK4")
xlabel("t [s]")
ylabel("T [K]")
legend("fixed step", "step doubling", "ode15s")
end
%% RK4 order and Test Equation
if mode == 4
%Find the order of RK4 by comparing the error of different fixed step sizes
%when approximating th etest equation up to t=1
global_error = [];
h = [];

% Get a nice amount of h sizes
for i=10:1:20
    [X, T] = RK4FixedStepSize([1.0], @testf, 1.0/i, 0, 1, 1);
    global_error = [global_error; abs(X(end)-exp(1))];
    h = [h;1.0/i];
end
for i=20:5:100
    [X, T] = RK4FixedStepSize([1.0], @testf, 1.0/i, 0, 1, 1);
    global_error = [global_error; abs(X(end)-exp(1))];
    h = [h;1.0/i];
end
for i=100:10:200
    [X, T] = RK4FixedStepSize([1.0], @testf, 1.0/i, 0, 1, 1);
    global_error = [global_error; abs(X(end)-exp(1))];
    h = [h;1.0/i];
end
for i=200:100:700
    [X, T] = RK4FixedStepSize([1.0], @testf, 1.0/i, 0, 1, 1);
    global_error = [global_error; abs(X(end)-exp(1))];
    h = [h;1.0/i];
end

%Compare the error to h^4, since this will show linearity if it is O(h^4)
plot(h.^4, global_error, 'Color', 'blue')
title("Global error vs h^4 for classical Runge Kutta")
xlabel("h^4")
ylabel("global error")
legend(["Global error"])
end
%% RK4 stability
if mode == 5
%Calculate |R(z)| in the square -4,-4 to 4,4. Everywhere |R(z)| < 1,
%absolute stability
x = -4:0.01:4;
y = -4:0.01:4;
v = @(z) abs(1+ z + 1/2*z^2 + 1/6*z^3 + 1/24*z^4);

value = zeros(length(y), length(x));
for n=1:length(x)
    for k=1:length(y)
        value(k,n) = v((x(n)+y(k)*1i));
    end
end
value(value>1) = 1;
image([-4,4], [-4,4], value, 'CDataMapping', 'scaled')
colorbar
title("|R(z)|, Classical Runge Kutta")
xlabel("Re(z)")
ylabel("Im(z)")

end