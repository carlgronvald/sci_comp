
global_error = [];
h = [];

for i=10:1:20
    [X, T] = RK4FixedStepSize([1.0], @testf, i, 0, 1, 1);
    global_error = [global_error; abs(X(end)-exp(1))];
    h = [h;1.0/i];
end
for i=20:5:100
    [X, T] = RK4FixedStepSize([1.0], @testf, i, 0, 1, 1);
    global_error = [global_error; abs(X(end)-exp(1))];
    h = [h;1.0/i];
end
for i=100:10:200
    [X, T] = RK4FixedStepSize([1.0], @testf, i, 0, 1, 1);
    global_error = [global_error; abs(X(end)-exp(1))];
    h = [h;1.0/i];
end
for i=200:100:700
    [X, T] = RK4FixedStepSize([1.0], @testf, i, 0, 1, 1);
    global_error = [global_error; abs(X(end)-exp(1))];
    h = [h;1.0/i];
end

hbase = global_error(end)/ h(end)^4;
plot(h.^4, global_error, 'Color', 'blue')
hold on
%plot(h, h.^4*0.019, 'Color', 'green')
hold off
title("Global error vs h^4 for classical Runge Kutta")
xlabel("h^4")
ylabel("global error")
legend(["Global error"])
saveas(gcf, "classical_rk_global_error.pdf")