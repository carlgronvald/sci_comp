function [Xres, fres] = NewtonsMethodESDIRK(X, T, f, h, gamma, psi, L, U, P, max_iter, epsilon, parameters)
%NEWTONSMETHODESDIRK Summary of this function goes here
%   Detailed explanation goes here
    R = ones(size(X,2), 1)
    X=X'
    iter = 0;
    while max(abs(R)) > epsilon
        fres = f(T, X, parameters);
        R = X - h*gamma*fres-psi;
        dx = -P*(U\(L\(P' * R)));
        X = X + dx;
        disp(R)
        
        iter = iter + 1;
        if iter > max_iter
            break;
        end
    end
    Xres = X
end

