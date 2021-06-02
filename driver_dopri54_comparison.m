% 1 == Van der Pol
mode = 3;
%% Van der Pol
if mode == 1
parameters = CreateParams('mu', 1.5);
x0 = [1.0;1.0];
vanmu1p5 = @(t,x) vanderpolf(t,x,parameters);
[X1,T1] = Dopri54(x0, @vanderpolf, 0.5, 0, 40, 0.01, 0.01, parameters);
[Tcorrect, Xcorrect] = ode15s(vanmu1p5, [0 40], [1.0;1.0]);
hold off
plot(T1, X1(:,1))
hold on
plot(Tcorrect, Xcorrect(:,1))
title("Van der Pol, mu=1.5, DOPRI54")
xlabel("t")
ylabel("value")
legend("x(1), DOPRI54 adaptive step", "x(1), ode15s")
figure 

parameters = CreateParams('mu', 15);
x0 = [1.0;1.0];
vanmu1p5 = @(t,x) vanderpolf(t,x,parameters);
[X1,T1] = Dopri54(x0, @vanderpolf, 0.5, 0, 40, 0.01, 0.01, parameters);
[Tcorrect, Xcorrect] = ode15s(vanmu1p5, [0 40], [1.0;1.0]);
hold off
plot(T1, X1(:,1))
hold on
plot(T2, X2(:,1))
plot(Tcorrect, Xcorrect(:,1))
title("Van der Pol, mu=15, DOPRI54")
xlabel("t")
ylabel("value")
legend("x(1), DOPRI54", "x(1), ode15s")
end
%% Adiabatic CSTR 3D
if mode == 2
parameters = CSTRparameters();
x0 = CSTRx0(parameters);
parmcstr = @(t,x) CSTRf(t,x,parameters);
[X1,T1] = Dopri54(x0, @CSTRf, 10, 0, 200, 0.01, 0.01, parameters);
[Tcorrect, Xcorrect] = ode15s(parmcstr, [0 200], x0);

hold off
plot(T1, X1(:,3))
hold on
plot(Tcorrect, Xcorrect(:,3))
xlabel("t [s]")
ylabel("T [K]")
title("3D CSTR, temperature, DOPRI54")
xlabel("t [s]")
ylabel("T [K]")
legend("DOPRI54", "ode15s")

end
%% Adiabatic CSTR 1D
if mode == 3
parameters = CSTRparameters();
x0 = 273.15;
parmcstr = @(t,x) CSTR1Df(t,x,parameters);
[X1,T1] = Dopri54(x0, @CSTR1Df, 10, 0, 200, 0.01, 0.01, parameters);
[Tcorrect, Xcorrect] = ode15s(parmcstr, [0 200], x0);

hold off
plot(T1, X1)
hold on
plot(Tcorrect, Xcorrect)
xlabel("t [s]")
ylabel("T [K]")
title("1D CSTR, temperature, DOPRI54")
xlabel("t [s]")
ylabel("T [K]")
legend("DOPRI54", "ode15s")
end
