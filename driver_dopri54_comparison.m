% 1 == Van der Pol
mode = 1;
%% Van der Pol
if mode == 1
parameters = CreateParams('mu', 1.5);
x0 = [1.0;1.0];

% Test DOPRI54 vs ode45, compare to ode15s. 
% Counter is used to count how many times ode45 and dopri54 call f
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

%We need to do the Van der Pol for mu=15 as well
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