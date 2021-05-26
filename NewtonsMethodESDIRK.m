function [Xres, fres, divergent] = NewtonsMethodESDIRK(X, T, f, h, gamma, psi, L, U, P, max_iter, epsilon, parameters)
%NEWTONSMETHODESDIRK Summary of this function goes here
%   Detailed explanation goes here
    R = Inf(size(X,1), 1);
    iter = 0;
    divergent = false;
    Xres = X;
    while max(abs(R)) > epsilon
        fres = f(T, X, parameters);
        Rlast = R;
        R = X - h*gamma*fres-psi;
        alpha = max(abs(R))/max(abs(Rlast));
        if alpha>1
            divergent = true;
            return; 
        end
        
        dx = -P*(U\(L\(P' * R)));
        X = X + dx;
        
        iter = iter + 1;
        if iter > max_iter
            divergent = true; %Slow convergence is as bad as divergence
            return;
        end
    end
    Xres = X;
end

