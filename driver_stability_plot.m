%% Esdirk Stability
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

%% RK4 Stability
x = -4:0.01:4;
y = -4:0.01:4;
v = @(z) abs(1+ z + 1/2*z^2 + 1/6*z^3 + 1/24*z^4);

value = zeros(length(y), length(x));
for n=1:length(x)
    for k=1:length(y)
        value(k,n) = v((x(n)+y(k)*1i));
    end
end
d = 255/5;
value(value>1) = 1;
image([-4,4], [-4,4], value, 'CDataMapping', 'scaled')
colorbar
title("|R(z)|, Classical Runge Kutta")
xlabel("Re(z)")
ylabel("Im(z)")

%% DOPRI54 Stability
Butcher = [0 1/5 3/10 4/5 8/9 1 1 0 0; ... 
    0 1/5 3/40 44/45 19372/6561 9017/3168 35/384 5179/57600 71/57600; ... 
    0 0 9/40 -56/15 -25360/2187 -355/33 0 0 0; ...
    0 0 0 32/9 64448/6561 46732/5247 500/1113 7571/16695 -71/16695; ...
    0 0 0 0 -212/729 49/176 125/192 393/640 71/1920; ...
    0 0 0 0 0 -5103/18656 -2187/6784 -92097/339200 -17253/339200; ...
    0 0 0 0 0 0 11/84 187/2100 22/525; ...
    0 0 0 0 0 0 0 1/40 -1/40]';
A = Butcher(1:7,2:8);
bf = Butcher(8,2:8);
% To calculate forward integrator:
% syms z
% calculate abs(1+z*bf*((eye(7)-z*A)^(-1))*[1;1;1;1;1;1;1]);
Rz = @(z) abs(1.0000*z - 0.6141*z*(- 0.8000*z^2 + 3.7333*z) - 0.2715*z*(- 0.0465*z^4 + 0.5816*z^3 - 1.8668*z^2 + 2.9526*z) + 0.0250*z*(0.0370*z^4 + 0.1111*z^3 + 0.3145*z^2 + 0.4492*z) - 0.0250*z*(- 0.0083*z^5 + 0.0139*z^4) + 0.0250*z*(0.0104*z^3 + 0.1302*z^2 + 0.6510*z) + 0.0890*z*(0.2828*z^3 - 1.6970*z^2 + 8.9064*z) + 0.2715*z*(0.2326*z^3 - 3.2958*z^2 + 11.5958*z) + 0.0250*z*(0.0017*z^6 + 0.0185*z^4 + 0.0451*z^3 + 0.0911*z^2 + 0.0911*z) + 0.6141*z*(0.1600*z^3 - 0.4800*z^2 + 0.9778*z) + 0.0890*z*(0.0127*z^5 - 0.1145*z^4 + 0.7778*z^3 - 2.0189*z^2 + 2.8463*z) + 2.3432*z^2 - 0.2715*z*(- 1.0340*z^2 + 9.8229*z) + 0.0890*z*(0.0795*z^2 + 0.2784*z) - 0.0250*z*(0.0358*z^2 + 0.3224*z) - 0.0890*z*(- 0.0636*z^4 + 0.6788*z^3 - 4.1364*z^2 + 10.7576*z) + 0.4535*z*(0.0450*z^2 + 0.0750*z) + 1);

x = -6:0.01:6;
y = -6:0.01:6;

value = zeros(length(y), length(x));
for n=1:length(x)
    for k=1:length(y)
        value(k,n) = Rz((x(n)+y(k)*1i));
    end
end
d = 255/5;
value(value>1) = 1;
image([x(1),x(end)], [y(1),y(end)], value, 'CDataMapping', 'scaled')
colorbar
title("|R(z)|, DOPRI54 Forward Integrator")
xlabel("Re(z)")
ylabel("Im(z)")

% Embedded Error Estimator
figure
be = Butcher(9,2:8);
Rz2 = @(z) abs(7.1222e-19*z - 0.0370*z*(- 0.8000*z^2 + 3.7333*z) - 0.0509*z*(- 0.0465*z^4 + 0.5816*z^3 - 1.8668*z^2 + 2.9526*z) - 0.0250*z*(0.0370*z^4 + 0.1111*z^3 + 0.3145*z^2 + 0.4492*z) + 0.0250*z*(- 0.0083*z^5 + 0.0139*z^4) - 0.0250*z*(0.0104*z^3 + 0.1302*z^2 + 0.6510*z) + 0.0419*z*(0.2828*z^3 - 1.6970*z^2 + 8.9064*z) + 0.0509*z*(0.2326*z^3 - 3.2958*z^2 + 11.5958*z) - 0.0250*z*(0.0017*z^6 + 0.0185*z^4 + 0.0451*z^3 + 0.0911*z^2 + 0.0911*z) + 0.0370*z*(0.1600*z^3 - 0.4800*z^2 + 0.9778*z) + 0.0419*z*(0.0127*z^5 - 0.1145*z^4 + 0.7778*z^3 - 2.0189*z^2 + 2.8463*z) + 0.1306*z^2 - 0.0509*z*(- 1.0340*z^2 + 9.8229*z) + 0.0419*z*(0.0795*z^2 + 0.2784*z) + 0.0250*z*(0.0358*z^2 + 0.3224*z) - 0.0419*z*(- 0.0636*z^4 + 0.6788*z^3 - 4.1364*z^2 + 10.7576*z) - 0.0043*z*(0.0450*z^2 + 0.0750*z) + 1);

value = zeros(length(y), length(x));
for n=1:length(x)
    for k=1:length(y)
        value(k,n) = Rz2((x(n)+y(k)*1i));
    end
end
d = 255/5;
value(value>1) = 1;
image([x(1),x(end)], [y(1),y(end)], value, 'CDataMapping', 'scaled')
colorbar
title("|R(z)|, DOPRI54 Embedded Error")
xlabel("Re(z)")
ylabel("Im(z)")
