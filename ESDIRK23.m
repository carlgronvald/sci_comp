function [X,T] = ESDIRK23(x0, f, jac, t0, t1, h0, abstol, reltol, parameters)
%ESDIRK23 Does some partially implicit partially explicit stuff I think
%jac = function that evaluates dx, jacobian
%
    epsilon = 0.8;
    gamma = (2-sqrt(2))/2;
    max_iter = 100;
    num_variables = size(x0, 1);
    
    bs = [ (1-gamma)/2; ...
        (1-gamma)/2; ...
        gamma];
    bhats = [(6*gamma-1)/(12*gamma); ...
        1/(12*gamma*(1-2*gamma)); ...
        (1-3*gamma)/(3*(1-2*gamma))];
    cs = [0;2*gamma;1];
    ds = [(1-6*gamma^2)/(12*gamma); ...
        (6*gamma*(1-2*gamma)*(1-gamma)-1)/(12*gamma*(1-2*gamma)); ...
        (6*gamma*(1-gamma)-1)/(3*(1-2*gamma))];
    divcount = 0;
    t = t0;
    x = x0;
    X = x0;
    T = t0;
    h = h0;
    hlast = h0;
    rlast = 0;
    Ts = zeros(1, 3);
    Xs = zeros(size(x0,1), 3);
    fs = zeros(size(x0,1), 3);
    

    fs(:,3) = f(t, x, parameters);
    n = 1;
    while t < t1
        if t+h > t1
           h = t1-t;
        end
        J = jac(t, x, parameters);
        fs(:,1) = fs(:,3);
        Ts(1) = t;
        Xs(:,1) = x;
        M = eye(size(J, 1))-h*gamma*J;
        [L, U, P] = lu(M);
        divergent = true;
        
        for i=2:3
            psi = x+h*( fs(:,1:i-1) * bs(1:i-1));
            Ts(i) = t + cs(i)*h;
            Xs(:,i) = x + cs(i)*h*fs(:,1);
            divergent = true;
            while divergent
                [Xres, fres, divergent] = NewtonsMethodESDIRK(Xs(:,i), Ts(i), f, h, gamma, psi, L, U, P, 20, 0.00000001, parameters);
                if divergent
                    disp("Divergence or slow convergence! Refactorizing...");
                    disp(["H=", h, ", t=", t])
                    divcount = divcount + 1;
                    h = h/2;
                    break;
                end
            end
            if divergent
                break;
            end
            Xs(:,i) = Xres;
            fs(:,i) = fres;
        end
        
        if divergent
            continue;
        end
        
        e = h*fs*ds;
        r = max(abs(e) ./ (abstol + abs(x)*reltol));
        if r<=1
            tmp = h; %Save old value of h
            if n==1 %Asymptotic error control for first step
                h = (epsilon/r)^(1/3)*h;
            else %Otherwise, PI error control
                h = (h/hlast)*(epsilon/r)^(1/3)*(rlast/r)^(1/3)*h;
            end
            x = Xs(:,3);
            t = Ts(3);
            X = [X x];
            T = [T;t];
            n = n+1;
            rlast = r; hlast = tmp;
        else
            h = (epsilon/r)^(1/3)*h;
        end
    end
    X = X';
    T = T';
    disp("divcount")
    disp(divcount)
end

