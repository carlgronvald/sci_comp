%1 == predator-prey
%2 == Van der Pol explicit euler
%3 == CSTR3D explicit euler
%4 == CSTR1D explicit euler
%5 == Van der Pol implicit Euler
%6 == CSTR3D implicit Euler
%7 == CSTR1D implicit euler
mode = 5;

%% Predator-Prey and Stability
if mode == 1
parameters = CreateParams('a', 1.1, 'b', 0.3);
x0 = [2;2];
[X1,T1] = ExplicitEulerFixedStepSize(x0, @predpreyf, 0.2, 0, 80, parameters);
[X2,T2] = ImplicitEulerFixedStepSize(x0, @predpreyf, @predpreyjac, 0.2, 0, 80, parameters);

hold off
plot(T2,X2(:,1))
hold on
plot(T2,X2(:,2))
title("The predator prey problem, implicit Euler")
xlabel("t")
ylabel("amount")
legend("prey", "predators")
figure
hold off
plot(T1, X1(:,1))
hold on
plot(T1,X1(:,2))
title("The predator prey problem, explicit Euler")
xlabel("t")
ylabel("amount")
legend("prey", "predators")
end
%% Van der Pol
if mode==2
parameters = CreateParams('mu', 1.5);
x0 = [1.0;1.0];
vanmu1p5 = @(t,x) vanderpolf(t,x,parameters);
[X1,T1] = ExplicitEulerFixedStepSize(x0, @vanderpolf, 0.2, 0, 40, parameters);
[X2,T2] = ExplicitEulerStepDoubling(x0, @vanderpolf, 0.2, 0, 40, 0.01, 0.01, parameters);
[Tcorrect, Xcorrect] = ode15s(vanmu1p5, [0 40], [1.0;1.0]);
hold off
plot(T1, X1(:,1))
hold on
plot(T2, X2(:,1))
plot(Tcorrect, Xcorrect(:,1))
title("Van der Pol, mu=1.5, explicit Euler")
xlabel("t")
ylabel("value")
legend("x(1), EE fixed step", "x(1), EE adaptive step", "x(1), ode15s")
figure 

parameters = CreateParams('mu', 12);
x0 = [1.0;1.0];
vanmu1p5 = @(t,x) vanderpolf(t,x,parameters);
[X1,T1] = ExplicitEulerFixedStepSize(x0, @vanderpolf, 0.01, 0, 40, parameters);
[X2,T2] = ExplicitEulerStepDoubling(x0, @vanderpolf, 0.2, 0, 40, 0.1, 0.1, parameters);
[Tcorrect, Xcorrect] = ode15s(vanmu1p5, [0 40], [1.0;1.0]);
hold off
plot(T1, X1(:,1))
hold on
plot(T2, X2(:,1))
plot(Tcorrect, Xcorrect(:,1))
title("Van der Pol, mu=15, explicit Euler")
xlabel("t")
ylabel("value")
legend("x(1), EE fixed step", "x(1), EE adaptive step", "x(1), ode15s")

end
%% CSTR explicit euler, 3D
if mode==3
parameters = CSTRparameters();
x0 = CSTRx0(parameters);
parmcstr = @(t,x) CSTRf(t,x,parameters);
[X1,T1] = ExplicitEulerFixedStepSize(x0, @CSTRf, 10, 0, 200, parameters);
[X2,T2] = ExplicitEulerStepDoubling(x0, @CSTRf, 1, 0, 200, 0.01, 0.01, parameters);
[Tcorrect, Xcorrect] = ode15s(parmcstr, [0 200], x0);

hold off
plot(T1, X1(:,3))
hold on
plot(T2, X2(:,3))
plot(Tcorrect, Xcorrect(:,3))
xlabel("t [s]")
ylabel("T [K]")
title("3D CSTR, temperature, explicit Euler")
xlabel("t [s]")
ylabel("T [K]")
legend("fixed step", "step doubling", "ode15s")
end
%% CSTR1D explicit euler
if mode==4
parameters = CSTRparameters();
x0 = 273.15;
parmcstr = @(t,x) CSTR1Df(t,x,parameters);
[X1,T1] = ExplicitEulerFixedStepSize(x0, @CSTR1Df, 10, 0, 200, parameters);
[X2,T2] = ExplicitEulerStepDoubling(x0, @CSTR1Df, 1, 0, 200, 0.01, 0.01, parameters);
[Tcorrect, Xcorrect] = ode15s(parmcstr, [0 200], x0);

hold off
plot(T1, X1)
hold on
plot(T2, X2)
plot(Tcorrect, Xcorrect)
xlabel("t [s]")
ylabel("T [K]")
title("1D CSTR, temperature, explicit Euler")
xlabel("t [s]")
ylabel("T [K]")
legend("fixed step", "step doubling", "ode15s")
end
%% Van der Pol implicit euler
if mode == 5
parameters = CreateParams('mu', 2);
x0 = [0.5;0.5];
vanmu1p5 = @(t,x) vanderpolf(t,x,parameters);
[X1,T1] = ImplicitEulerFixedStepSize(x0, @vanderpolf, @vanderpoljac, 0.2, 0, 40, parameters);
[X2,T2] = ImplicitEulerStepDoubling(x0, @vanderpolf, @vanderpoljac, 0.2, 0, 40, 0.01, 0.01, parameters);
[Tcorrect, Xcorrect] = ode15s(vanmu1p5, [0 40], [1.0;1.0]);
hold off
plot(T1, X1(:,1))
hold on
plot(T2, X2(:,1))
plot(Tcorrect, Xcorrect(:,1))
title("Van der Pol, mu=2, implicit Euler")
xlabel("t")
ylabel("value")
legend("x(1), IE fixed step", "x(1), IE adaptive step", "x(1), ode15s")
figure 

parameters = CreateParams('mu', 12.0);
x0 = [1.0;1.0];
vanmu1p5 = @(t,x) vanderpolf(t,x,parameters);
[X1,T1] = ImplicitEulerFixedStepSize(x0, @vanderpolf, @vanderpoljac, 0.05, 0, 40, parameters);
[X2,T2] = ImplicitEulerStepDoubling(x0, @vanderpolf, @vanderpoljac, 0.2, 0, 40, 0.005, 0.005, parameters);
[Tcorrect, Xcorrect] = ode15s(vanmu1p5, [0 40], [1.0;1.0]);
hold off
plot(T1, X1(:,1))
hold on
plot(T2, X2(:,1))
plot(Tcorrect, Xcorrect(:,1))
title("Van der Pol, mu=12, implicit Euler")
xlabel("t")
ylabel("value")
legend("x(1), IE fixed step", "x(1), IE adaptive step", "x(1), ode15s")
    
end

%% CSTR implicit euler, 3D
if mode==6
parameters = CSTRparameters();
x0 = CSTRx0(parameters);
parmcstr = @(t,x) CSTRf(t,x,parameters);
[X1,T1] = ImplicitEulerFixedStepSize(x0, @CSTRf, @CSTRjac, 10, 0, 200, parameters);
[X2,T2] = ImplicitEulerStepDoubling(x0, @CSTRf, @CSTRjac, 1, 0, 200, 0.01, 0.01, parameters);
[Tcorrect, Xcorrect] = ode15s(parmcstr, [0 200], x0);

hold off
plot(T1, X1(:,3))
hold on
plot(T2, X2(:,3))
plot(Tcorrect, Xcorrect(:,3))
xlabel("t [s]")
ylabel("T [K]")
title("3D CSTR, temperature, implicit Euler")
xlabel("t [s]")
ylabel("T [K]")
legend("fixed step", "step doubling", "ode15s")
end
%% CSTR1D implicit euler
if mode==7
parameters = CSTRparameters();
x0 = 273.15;
parmcstr = @(t,x) CSTR1Df(t,x,parameters);
[X1,T1] = ImplicitEulerFixedStepSize(x0, @CSTR1Df, @CSTR1Djac, 10, 0, 200, parameters);
[X2,T2] = ImplicitEulerStepDoubling(x0, @CSTR1Df, @CSTR1Djac, 1, 0, 200, 0.01, 0.01, parameters);
[Tcorrect, Xcorrect] = ode15s(parmcstr, [0 200], x0);

hold off
plot(T1, X1)
hold on
plot(T2, X2)
plot(Tcorrect, Xcorrect)
xlabel("t [s]")
ylabel("T [K]")
title("1D CSTR, temperature, implicit Euler")
xlabel("t [s]")
ylabel("T [K]")
legend("fixed step", "step doubling", "ode15s")
end
