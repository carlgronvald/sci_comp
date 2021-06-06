% 1 == Van der Pol ESDIRK23
% 2 == CSTR 3D ESDIRK23
% 3 == CSTR 1D ESDIRK23
% 4 == Extreme Van der Pol, mu=500
% 5 == ESDIRK23 forward integrator stability
mode = 4;
%% Van der Pol
if mode == 1
parameters = CreateParams('mu', 1.5);
x0 = [1.0;1.0];
vanmu1p5 = @(t,x) vanderpolf(t,x,parameters);
[X1,T1] = ESDIRK23(x0, @vanderpolf, @vanderpoljac, 0.01, 0, 40, 0.01, 0.01, parameters);
[Tcorrect, Xcorrect] = ode15s(vanmu1p5, [0 40], [1.0;1.0]);
hold off
plot(T1, X1(:,1))
hold on
plot(Tcorrect, Xcorrect(:,1))
title("Van der Pol, mu=1.5, ESDIRK23")
xlabel("t")
ylabel("x(1)")
legend("ESDIRK23 adaptive step", "ode15s")
figure 

parameters = CreateParams('mu', 15);
x0 = [1.0;1.0];
vanmu1p5 = @(t,x) vanderpolf(t,x,parameters);
[X1,T1] = ESDIRK23(x0, @vanderpolf, @vanderpoljac, 0.01, 0, 40, 0.01, 0.01, parameters);
[Tcorrect, Xcorrect] = ode15s(vanmu1p5, [0 40], [1.0;1.0]);
hold off
plot(T1, X1(:,1))
hold on
plot(Tcorrect, Xcorrect(:,1))
title("Van der Pol, mu=15, ESDIRK23")
xlabel("t")
ylabel("x(1)")
legend("ESDIRK23", "ode15s")
end
%% Adiabatic CSTR 3D
if mode == 2
parameters = CSTRparameters();
x0 = CSTRx0(parameters);
parmcstr = @(t,x) CSTRf(t,x,parameters);
[X1,T1] = ESDIRK23(x0, @CSTRf, @CSTRjac, 0.01, 0, 200, 0.01, 0.01, parameters);
[Tcorrect, Xcorrect] = ode15s(parmcstr, [0 200], x0);

hold off
plot(T1, X1(:,3))
hold on
plot(Tcorrect, Xcorrect(:,3))
xlabel("t [s]")
ylabel("T [K]")
title("3D CSTR, temperature, ESDIRK23")
xlabel("t [s]")
ylabel("T [K]")
legend("ESDIRK23", "ode15s")

end
%% Adiabatic CSTR 1D
if mode == 3
parameters = CSTRparameters();
x0 = 273.15;
parmcstr = @(t,x) CSTR1Df(t,x,parameters);
[X1,T1] = ESDIRK23(x0, @CSTR1Df, @CSTR1Djac, 0.01, 0, 200, 0.001, 0.001, parameters);
[Tcorrect, Xcorrect] = ode15s(parmcstr, [0 200], x0);

hold off
plot(T1, X1)
hold on
plot(Tcorrect, Xcorrect)
xlabel("t [s]")
ylabel("T [K]")
title("1D CSTR, temperature, ESDIRK23")
xlabel("t [s]")
ylabel("T [K]")
legend("ESDIRK23", "ode15s")
end
%% Extreme Van der Pol
if mode == 4
parameters = CreateParams('mu', 500);
x0 = [1.0;1.0];
vanmu1p5 = @(t,x) vpcounter(t,x,parameters);
% We want to measure both the number of function evaluations and the number
% of jacobian evaluations. Only ESDIRK23 uses jacobians, the others only
% use function evaluations.
% We also measure the time used by each method.
global counter;
global Jcounter;
counter=0;
Jcounter=0;

start = cputime;
[X1,T1] = ESDIRK23(x0, @vpcounter, @vpjcounter, 0.01, 0, 800, 0.01, 0.01, parameters);
disp(cputime-start);
disp(["counter", counter, "JCounter", Jcounter]);

counter =0;
start = cputime;
[X2,T2] = Dopri54(x0, @vpcounter, 0.01, 0, 800, 0.01, 0.01, parameters);
disp(cputime-start);
disp(["counter", counter])

counter = 0;
start = cputime;
[X3,T3] = RK4StepDoubling(x0, @vpcounter, 0.01, 0, 800, 0.01, 0.01, parameters);
disp(cputime-start)
disp(["counter", counter])
counter=0;
tic
[Tcorrect, Xcorrect] = ode15s(vanmu1p5, [0 800], [1.0;1.0]);
disp(toc)
disp(["counter", counter])

hold off
plot(T1, X1(:,1))
hold on
plot(T2, X2(:,1))
plot(T3, X3(:,1))
plot(Tcorrect, Xcorrect(:,1))
title("Van der Pol, mu=500")
xlabel("t")
ylabel("x(1)")
legend("ESDIRK23", "DOPRI54", "RK4 step doubling", "ode15s")
end

%% Esdirk Stability
if mode == 5
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
end


function dx = vpcounter(t,x,p)
    global counter;
    dx= vanderpolf(t,x,p);
    counter = counter+1;
end

function J = vpjcounter(t,x,p)
    global Jcounter;
    J = vanderpoljac(t,x,p);
    Jcounter = Jcounter+1;
end    
