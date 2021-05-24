function [outputArg1,outputArg2] = ESDIRK23(x0, f, jac, t0, t1, h0, parameters)
%ESDIRK23 Does some partially implicit partially explicit stuff I think
%jac = function that evaluates dx, jacobian
%
    gamma = (2-sqrt(2))/2;
    max_iter = 100;
    
    bs = [(1-gamma)/2;(1-gamma)/2;gamma];
    bhats = [(6*gamma-1)/(12*gamma); ...
        1/(12*gamma*(1-2*gamma)); ...
        (1-3*gamma)/(3*(1-2*gamma))];
    cs = [0;2*gamma;1];
    t = t0;
    x = x0;
    h = h0;
    Ts = zeros(3, 1);
    Xs = zeros(3, size(x0,1));
    fs = zeros(3, size(x0,1));

    fs(3,:) = f(t, x);
    while t < t1
        
        J = jac(t, x, parameters);
        fs(1,:) = fs(3,:);
        Ts(1) = t;
        Xs(1,:) = x;
        M = eye(size(J, 1))-h*gamma*J;
        [L, U, P] = lu(M);
        
        for i=2:3
            psi = x+h*( reshape(bs(1:i-1),[1,i-1]) * fs(1:i-1,:));
            Ts(i) = t + cs(i)*h;
            Xs(i,:) = x + cs(i)*h*fs(1);
            [Xres, fres] = NewtonsMethodESDIRK(Xs(i,:), Ts(i), f, h, gamma, psi, L, U, P, 100, 0.001, parameters);
            Xs(i,:) = Xres;
            fs(i,:) = fres;
        end
        x = Xs(3,:)
        %TODO: OPTIMIZATIONS
    end
end

