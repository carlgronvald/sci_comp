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