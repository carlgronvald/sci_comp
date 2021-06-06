mode = 1;

if mode == 1
W = Wiener(0, 40, 200, 2, 1, 1);
p = CreateParams('mu', 1.5, 'sigma', 0.3);
[X1, T1] = SDEExplicitExplicitFixedStepSize([1.0;1.0], @vanderpolf, @vanderpolg, 0.2, 0, 40, W, p);
[X2, T2] = SDEImplicitExplicitFixedStepSize([1.0;1.0], @vanderpolf, @vanderpoljac, @vanderpolg, 0.2, 0, 40, W, p);

hold off
plot(T1, X1(:,1));
hold on
plot(T2, X2(:,1));
xlabel("t")
ylabel("x(1)")
legend("Explicit-Explicit", "Implicit-Explicit")
title("Stochastic Van der Pol, mu=1.5, sigma=0.3")
figure

W = Wiener(0, 40, 4000, 2, 1, 1);
p = CreateParams('mu', 15, 'sigma', 1.5);
[X1, T1] = SDEExplicitExplicitFixedStepSize([1.0;1.0], @vanderpolf, @vanderpolg, 0.01, 0, 40, W, p);
[X2, T2] = SDEImplicitExplicitFixedStepSize([1.0;1.0], @vanderpolf, @vanderpoljac, @vanderpolg, 0.01, 0, 40, W, p);

hold off
plot(T1, X1(:,1));
hold on
plot(T2, X2(:,1));
xlabel("t")
ylabel("x(1)")
legend("Explicit-Explicit", "Implicit-Explicit")
title("Stochastic Van der Pol, mu=15, sigma=1.5")
end
%% CSTR 3D
if mode == 2
W = Wiener(0, 200, 200, 3, 1, 1);
p = CSTRparameters('sigma', 1.5);
x0 = CSTRx0(p);
[X1, T1] = SDEExplicitExplicitFixedStepSize(x0, @CSTRf, @CSTRg, 1, 0, 200, W, p);
[X2, T2] = SDEImplicitExplicitFixedStepSize(x0, @CSTRf, @CSTRjac, @CSTRg, 1, 0, 200, W, p);

plot(T1, X1(:,3))
hold on
plot(T2, X2(:,3))
hold off
xlabel("t [s]")
ylabel("T [K]")
legend("Explicit-Explicit", "Implicit-Explicit")
title("Stochastic 3D CSTR, sigma=1.5")
end
%% CSTR 1D
if mode == 3
W = Wiener(0, 200, 200, 1, 1, 1);
p = CSTRparameters('sigma', 1.5);
x0 = [273.15];
[X1, T1] = SDEExplicitExplicitFixedStepSize(x0, @CSTR1Df, @CSTR1Dg, 1, 0, 200, W, p);
[X2, T2] = SDEImplicitExplicitFixedStepSize(x0, @CSTR1Df, @CSTR1Djac, @CSTR1Dg, 1, 0, 200, W, p);

plot(T1, X1)
hold on
plot(T2, X2)
hold off
xlabel("t [s]")
ylabel("T [K]")
title("Stochastic 1D CSTR, sigma=1.5")
legend("Explicit-Explicit", "Implicit-Explicit")
end