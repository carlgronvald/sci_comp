% 1 == Van der Pol DOPRI54
% 2 == 3D CSTR DOPRI54
% 3 == 1D CSTR DOPRI54
% 4 == DOPRI54 Stability
mode = 4;
%% Van der Pol
if mode == 1
parameters = CreateParams('mu', 1.5);
x0 = [1.0;1.0];

% Test DOPRI54 vs ode45, compare to ode15s. 
% Counter is used to count how many times ode45 and dopri54 call f
% Then, we can compare how efficient our solution is compared to the Matlab
% solution.
global counter;
counter = 0;
vanmu1p5 = @(t,x) vpcounter(t,x,parameters);
[X1,T1] = Dopri54(x0, @vpcounter, 1, 0, 40, 0.01, 0.01, parameters);
disp(counter)
counter = 0;
%ode45 and dopri54 are given same relative and absolute tolerance
options = odeset('RelTol', 0.01, 'AbsTol', 0.01);
[T2,X2] = ode45(vanmu1p5, [0 40], x0, options);
disp(counter)
[Tcorrect, Xcorrect] = ode15s(vanmu1p5, [0 40], x0);

%Plot them together.
hold off
plot(T1, X1(:,1))
hold on
plot(T2, X2(:,1))
plot(Tcorrect, Xcorrect(:,1))
title("Van der Pol, mu=1.5, DOPRI54")
xlabel("t")
ylabel("x(1)")
legend("DOPRI54", "ode45", "ode15s")
figure 

%Now, Van der Pol with mu=15
parameters = CreateParams('mu', 15);
x0 = [1.0;1.0];
vanmu1p5 = @(t,x) vpcounter(t,x,parameters);
counter = 0;
[X1,T1] = Dopri54(x0, @vpcounter, 0.5, 0, 40, 0.01, 0.01, parameters);
disp(counter);
counter=0;
%ode45 and dopri54 are given same relative and absolute tolerance
options = odeset('RelTol', 0.01, 'AbsTol', 0.01);
[T2,X2] = ode45(vanmu1p5, [0 40], x0, options);
disp(counter);
[Tcorrect, Xcorrect] = ode15s(vanmu1p5, [0 40], [1.0;1.0]);

%Plot them together
hold off
plot(T1, X1(:,1))
hold on
plot(T2, X2(:,1))
plot(Tcorrect, Xcorrect(:,1))
title("Van der Pol, mu=15, DOPRI54")
xlabel("t")
ylabel("x(1)")
legend("DOPRI54", "ode45", "ode15s")

end
%% Adiabatic CSTR 3D
if mode == 2
%We solve the Adiabatic CSTR 3D problem using DOPRI54
%Again, we compare the number of function evaluations of DOPRI54 and ode45
parameters = CSTRparameters();
x0 = CSTRx0(parameters);
parmcstr = @(t,x) c3counter(t,x,parameters);
global counter;
counter = 0;
[X1,T1] = Dopri54(x0, @c3counter, 0.01, 0, 200, 0.01, 0.01, parameters);
disp(counter);
counter=0;
options = odeset('RelTol', 0.01, 'AbsTol', 0.01);
[T2,X2] = ode45(parmcstr, [0 200], x0, options);
disp(counter)
[Tcorrect, Xcorrect] = ode15s(parmcstr, [0 200], x0);

hold off
plot(T1, X1(:,3))
hold on
plot(T2, X2(:,3))
plot(Tcorrect, Xcorrect(:,3))
xlabel("t [s]")
ylabel("T [K]")
title("3D CSTR, temperature, DOPRI54")
xlabel("t [s]")
ylabel("T [K]")
legend("DOPRI54", "ode45", "ode15s")

end
%% Adiabatic CSTR 1D
if mode == 3
%We solve CSTR 1D using Dopri54
parameters = CSTRparameters();
x0 = 273.15;

% Again, we count the number of function evaluations
% We compare to ode15s and ode45
global counter;
counter = 0;
parmcstr = @(t,x) c1counter(t,x,parameters);
[X1,T1] = Dopri54(x0, @c1counter, 0.01, 0, 200, 0.01, 0.01, parameters);
disp(counter)
counter = 0;
options = odeset('RelTol', 0.01, 'AbsTol', 0.01);
[T2,X2] = ode45(parmcstr, [0 200], x0, options);
disp(counter)
[Tcorrect, Xcorrect] = ode15s(parmcstr, [0 200], x0);


hold off
plot(T1, X1)
hold on
plot(T2, X2)
plot(Tcorrect, Xcorrect)
xlabel("t [s]")
ylabel("T [K]")
title("1D CSTR, temperature, DOPRI54")
legend("DOPRI54", "ode45", "ode15s")
end
%% DOPRI54 Stability
if mode == 4
%The stability region of DOPRI54 both forward integrator and embedded
%method is calculated using the Butcher tableau
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
% To calculate forward integrator stability:
% syms z
% calculate abs(1+z*bf*((eye(7)-z*A)^(-1))*[1;1;1;1;1;1;1])
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
%again, calculated by
% syms z
% abs(1+z*be*((eye(7)-z*A)^(-1))*[1;1;1;1;1;1;1])
Rz2 = @(z) abs(7.1222e-19*z - 0.0370*z*(- 0.8000*z^2 + 3.7333*z) - 0.0509*z*(- 0.0465*z^4 + 0.5816*z^3 - 1.8668*z^2 + 2.9526*z) - 0.0250*z*(0.0370*z^4 + 0.1111*z^3 + 0.3145*z^2 + 0.4492*z) + 0.0250*z*(- 0.0083*z^5 + 0.0139*z^4) - 0.0250*z*(0.0104*z^3 + 0.1302*z^2 + 0.6510*z) + 0.0419*z*(0.2828*z^3 - 1.6970*z^2 + 8.9064*z) + 0.0509*z*(0.2326*z^3 - 3.2958*z^2 + 11.5958*z) - 0.0250*z*(0.0017*z^6 + 0.0185*z^4 + 0.0451*z^3 + 0.0911*z^2 + 0.0911*z) + 0.0370*z*(0.1600*z^3 - 0.4800*z^2 + 0.9778*z) + 0.0419*z*(0.0127*z^5 - 0.1145*z^4 + 0.7778*z^3 - 2.0189*z^2 + 2.8463*z) + 0.1306*z^2 - 0.0509*z*(- 1.0340*z^2 + 9.8229*z) + 0.0419*z*(0.0795*z^2 + 0.2784*z) + 0.0250*z*(0.0358*z^2 + 0.3224*z) - 0.0419*z*(- 0.0636*z^4 + 0.6788*z^3 - 4.1364*z^2 + 10.7576*z) - 0.0043*z*(0.0450*z^2 + 0.0750*z) + 1);

value = zeros(length(y), length(x));
for n=1:length(x)
    for k=1:length(y)
        value(k,n) = Rz2((x(n)+y(k)*1i));
    end
end
value(value>1) = 1;
image([x(1),x(end)], [y(1),y(end)], value, 'CDataMapping', 'scaled')
colorbar
title("|R(z)|, DOPRI54 Embedded Error")
xlabel("Re(z)")
ylabel("Im(z)")
end


%these functions are counting functions to check how many function calls
%are made.
function dx = vpcounter(t,x,p)
    global counter;
    dx = vanderpolf(t,x,p);
    counter = counter+1;
end

function dx = c1counter(t,x,p)
    global counter;
    dx = CSTR1Df(t,x,p);
    counter = counter+1;
end

function dx = c3counter(t,x,p)
    global counter;
    dx = CSTRf(t,x,p);
    counter = counter+1;
end