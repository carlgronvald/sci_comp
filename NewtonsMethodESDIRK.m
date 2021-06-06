function [Xres, fres, divergent] = NewtonsMethodESDIRK(X, T, f, h, gamma, psi, L, U, P, max_iter, epsilon, parameters)
%NEWTONSMETHODESDIRK A modified Newton's method using some improvements for
%ESDIRK23. Exits early in case of divergence, and counts slow convergence
%as divergence
% Note that this method does not calculate the Jacobian each time, but
% instead uses the factorization of the iteration matrix M. This is one
% reason we would rather just retry if we do not converge well, since the
% issue is likely with the M matrix that needs refactorization
%
% X - Initial guess for this stage of ESDIRK23
% T - time step of this stage of ESDIRK23
% f - the ODE function
% h - The timestep
% gamma - gamma for the ESDIRK method
% psi - Remainder term in the expression for the resdiual
% L - L from the LU factorization of M
% U - U from the LU factorization of M
% P - permutation matrix from the LU factorization of M
% max_iter - maximum number of iterations
% epsilon - tolerance for size of the residual
% parameters - parameters for the ODE function
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

