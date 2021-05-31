
%% Van der Pol
parameters = CreateParams('mu', 1.5);
x0 = [1.0;1.0];
[X1,T1] = ExplicitEulerFixedStepSize(x0, @vanderpolf, 100, 0, 10, parameters);
[X2,T2] = ExplicitEulerStepDoubling(x0, @vanderpolf, 0.1, 0, 10, 0.01, 0.01, parameters);
[X3,T3] = ImplicitEulerFixedStepSize(x0, @vanderpolf, @vanderpoljac, 100, 0, 10, parameters);
[X4,T4] = ImplicitEulerStepDoubling(x0, @vanderpolf, @vanderpoljac, 0.1, 0, 10, 0.01, 0.01, parameters);
[X5,T5] = RK4FixedStepSize(x0, @vanderpolf, 100, 0, 10, parameters);
[X6,T6] = RK4StepDoubling(x0, @vanderpolf, 0, 10, 0.1, 0.01, 0.01, parameters);
[X7,T7] = Dopri54AdaptiveStepSize(x0, @vanderpolf, 0.1, 0, 10, 0.01, 0.01, parameters);
[X8,T8] = ESDIRK23(x0, @vanderpolf, @vanderpoljac, 0, 10, 0.1, 0.01, 0.01, parameters);
[Xcorrect,Tcorrect] = ESDIRK23(x0,@vanderpolf,@vanderpoljac,0,10,0.01,1e-7,1e-7,parameters);


gl1 = find_global_error(X1,T1,Xcorrect,Tcorrect);
gl2 = find_global_error(X2,T2,Xcorrect,Tcorrect);
gl3 = find_global_error(X3,T3,Xcorrect,Tcorrect);
gl4 = find_global_error(X4,T4,Xcorrect,Tcorrect);
gl5 = find_global_error(X5,T5,Xcorrect,Tcorrect);
gl6 = find_global_error(X6,T6,Xcorrect,Tcorrect);
gl7 = find_global_error(X7,T7,Xcorrect,Tcorrect);
gl8 = find_global_error(X8,T8,Xcorrect,Tcorrect);

hold off
%plot(T1,gl1(:,1))
%hold on
%plot(T2,gl2(:,1))
%plot(T3,gl3(:,1))
%plot(T4,gl4(:,1))
%plot(T5,gl5(:,1))
%plot(T6,gl6(:,1))
plot(T7,gl7(:,1))
hold on
plot(T8,gl8(:,1))
legend(["EE, fixed step", "EE, adaptive step", "IE, fixed step", "IE, adaptive step", "RK4C, fixed step","RK4C, adaptive step", "DOPRI5(4) adaptive step", "ESDIRK23 adaptive step"])
title("Error at time t")
xlabel("t")
ylabel("Absolute error of x(1)")

%hold off
%plot(X1(:,1),X1(:,2));
%hold on
%plot(X2(:,1),X2(:,2));
%plot(X3(:,1),X3(:,2));
%plot(X4(:,1),X4(:,2));
%plot(X5(:,1),X5(:,2));
%plot(X6(:,1),X6(:,2));
%plot(X7(:,1),X7(:,2));
%plot(X8(:,1),X8(:,2));
%plot(Xcorrect(:,1),Xcorrect(:,2));
legend(["EE, fixed step", "EE, adaptive step", "IE, fixed step", "IE, adaptive step", "RK4C, fixed step","RK4C, adaptive step", "DOPRI5(4) adaptive step", "ESDIRK23 adaptive step"])
title("Many approximations of the ODE Van der Pol problem")
xlabel("x1")
ylabel("x2")

%% Helper functions
function [X] = closest_step(T, Tcorrect, Xcorrect)
    d = Tcorrect-T;
    [minValue, closestIndex] = min(abs(d(d<0)));
    Tmin = Tcorrect(closestIndex);
    Tmax = Tcorrect(closestIndex+1);
    Tdiff = T-Tmin;
    alpha = Tdiff/(Tmax-Tmin);
    X = Xcorrect(closestIndex,:) * (1-alpha) + Xcorrect(closestIndex+1,:)*alpha;
    
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