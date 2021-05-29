
%% Van der Pol
parameters = CreateParams('mu', 1.5);
x0 = [1.0;1.0];
[X1,T1] = ExplicitEulerFixedStepSize(x0, @vanderpolf, 100, 0, 10, parameters);
[X2,T2] = ExplicitEulerStepDoubling(x0, @vanderpolf, 0.1, 0, 10, 0.01, 0.01, parameters);
[X3,T3] = ImplicitEulerFixedStepSize(x0, @vanderpolf, @vanderpoljac, 100, 0, 10, parameters);
[X4,T4] = ImplicitEulerStepDoubling(x0, @vanderpolf, @vanderpoljac, 0.1, 0, 10, 0.01, 0.01, parameters);
[X5,T5] = Dopri54AdaptiveStepSize(x0, @vanderpolf, 0.1, 0, 10, 0.01, 0.01, parameters)
%[X5,T5] = ESDIRK23(x0, @vanderpolf, @vanderpoljac, 0, 10, 0.1, parameters);

hold off
plot(X1(:,1),X1(:,2));
hold on
plot(X2(:,1),X2(:,2));
plot(X3(:,1),X3(:,2));
plot(X4(:,1),X4(:,2));
legend(["EE, fixed step", "EE, adaptive step", "IE, fixed step", "IE, adaptive step"])
title("Many approximations of the Van der Pol problem")
xlabel("x1")
ylabel("x2")