
%% Van der Pol

parameters = CreateParams('mu', 1.5);
x0 = [1.0;1.0];
[X1,T1] = ExplicitEulerFixedStepSize(x0, @vanderpolf, 0.0007, 0, 40, parameters);
[X2,T2] = ExplicitEulerStepDoubling(x0, @vanderpolf, 0.1, 0, 40, 0.0000035, 0.000001, parameters);
[X3,T3] = ImplicitEulerFixedStepSize(x0, @vanderpolf, @vanderpoljac, 0.0007, 0, 40, parameters);
[X4,T4] = ImplicitEulerStepDoubling(x0, @vanderpolf, @vanderpoljac, 0.1, 0, 40, 0.000003, 0.000003, parameters);
[X5,T5] = RK4FixedStepSize(x0, @vanderpolf, 0.24, 0, 40, parameters);
[X6,T6] = RK4StepDoubling(x0, @vanderpolf, 0.1, 0, 40, 0.198, 0.198, parameters);
[X7,T7] = Dopri54(x0, @vanderpolf, 0.1, 0, 40, 0.001001, 0.001001, parameters);
[X8,T8] = ESDIRK23(x0, @vanderpolf, @vanderpoljac, 0.1, 0, 40, 0.000002, 0.000002, parameters);
options = odeset('RelTol', 1e-15, 'AbsTol', 1e-15);
vanmu1p5 = @(t,x) vanderpolf(t,x,parameters);   
[Tcorrect, Xcorrect] = ode15s(vanmu1p5, [0 40], [1.0;1.0], options);

gl1 = max(find_global_error(X1, T1, Xcorrect, Tcorrect));
gl1 = gl1(1);
gl2 = max(find_global_error(X2, T2, Xcorrect, Tcorrect));
gl2 = gl2(1);
gl3 = max(find_global_error(X3, T3, Xcorrect, Tcorrect));
gl3 = gl3(1);
gl4 = max(find_global_error(X4, T4, Xcorrect, Tcorrect));
gl4 = gl4(1);
gl5 = max(find_global_error(X5, T5, Xcorrect, Tcorrect));
gl5 = gl5(1);
gl6 = max(find_global_error(X6, T6, Xcorrect, Tcorrect));
gl6 = gl6(1);
gl7 = max(find_global_error(X7, T7, Xcorrect, Tcorrect));
gl7 = gl7(1);
gl8 = max(find_global_error(X8, T8, Xcorrect, Tcorrect));
gl8 = gl8(1);


%gl1 = Xcorrect(end)-X1(end);
%gl2 = Xcorrect(end)-X2(end);
%gl3 = Xcorrect(end)-X3(end);
%gl4 = Xcorrect(end)-X4(end);
%gl5 = Xcorrect(end)-X5(end);
%gl6 = Xcorrect(end)-X6(end);
%gl7 = Xcorrect(end)-X7(end);
%gl8 = Xcorrect(end)-X8(end);
%gl1 = find_global_error(X1,T1,Xcorrect,Tcorrect);
%gl2 = find_global_error(X2,T2,Xcorrect,Tcorrect);
%gl3 = find_global_error(X3,T3,Xcorrect,Tcorrect);
%gl4 = find_global_error(X4,T4,Xcorrect,Tcorrect);
%gl5 = find_global_error(X5,T5,Xcorrect,Tcorrect);
%gl6 = find_global_error(X6,T6,Xcorrect,Tcorrect);
%gl7 = find_global_error(X7,T7,Xcorrect,Tcorrect);
%gl8 = find_global_error(X8,T8,Xcorrect,Tcorrect);

%gl1 = max(gl1(:,1));
%gl2 = max(gl2(:,1));
%gl3 = max(gl3(:,1));
%gl4 = max(gl4(:,1));
%gl5 = max(gl5(:,1));
%gl6 = max(gl6(:,1));
%gl7 = max(gl7(:,1));
%gl8 = max(gl8(:,1));


hold off
%plot(T1,gl1(:,1))
%plot(T2,gl2(:,1))
plot(T3,X3(:,1))
%plot(T4,gl4(:,1))
%plot(T5,gl5(:,1))
%plot(T6,gl6(:,1))
hold on
plot(T7,X7(:,1))
hold on
plot(T8,X8(:,1))
plot(Tcorrect,Xcorrect(:,1))
legend(["EE, fixed step", "EE, adaptive step", "IE, fixed step", "IE, adaptive step", "RK4C, fixed step","RK4C, adaptive step", "DOPRI5(4) adaptive step", "ESDIRK23 adaptive step", "Very exact ode15s"])
title("Error at time t")
xlabel("t")
ylabel("Absolute error of x(1)")

hold off
plot(X1(:,1),X1(:,2));
hold on
plot(X2(:,1),X2(:,2));
plot(X3(:,1),X3(:,2));
plot(X4(:,1),X4(:,2));
plot(X5(:,1),X5(:,2));
plot(X6(:,1),X6(:,2));
plot(X7(:,1),X7(:,2));
plot(X8(:,1),X8(:,2));
plot(Xcorrect(:,1),Xcorrect(:,2));
legend(["EE, fixed step", "EE, adaptive step", "IE, fixed step", "IE, adaptive step", "RK4C, fixed step","RK4C, adaptive step", "DOPRI5(4) adaptive step", "ESDIRK23 adaptive step"])
title("Many approximations of the ODE Van der Pol problem")
xlabel("x1")
ylabel("x2")

%% Helper functions
function [X] = closest_step(T, Tcorrect, Xcorrect)
    d = Tcorrect-T;
    [minValue, closestIndex] = min(abs(d));
%    disp("lol")
%    disp(closestIndex)
%    alpha = Tdiff/(Tmax-Tmin);
%    X = Xcorrect(closestIndex,:) * (1-alpha) + Xcorrect(closestIndex+1,:)*alpha;
    X = Xcorrect(closestIndex, :);
    
end

function global_error = find_global_error(X, T, Xcorrect, Tcorrect)

global_error = zeros(length(T), 2);
for i=1:length(T)
    Xc = closest_step(T(i), Tcorrect, Xcorrect);
    e=abs(X(i,:)-Xc);
    if length(e)>0
        global_error(i,:) = e;
    end
end

end