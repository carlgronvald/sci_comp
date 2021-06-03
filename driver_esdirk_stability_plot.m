x = -10:0.01:15;
y = -10:0.01:10;
g = (2 - sqrt(2))/2;
v = @(x,y) abs((1 + (1-2*g)*(x + y*1i))/( (1 - g * (x+y*1i))^2));

value = zeros(length(y), length(x));
for n=1:length(x)
    for k=1:length(y)
        value(k,n) = v(x(n),y(k));
    end
end
d = 255/5;
value(value>1) = 1;
image([-10,15], [-10,10], value, 'CDataMapping', 'scaled')
colorbar
title("|R(z)|, ESDIRK23 Forward Integrator")
xlabel("Re(z)")
ylabel("Im(z)")