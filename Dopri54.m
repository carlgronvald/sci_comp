function [X,T] = Dopri54(x0, f, h0, t0, t1, abstol, reltol, params)
%DOPRI54 Solves an ODE using the DOPRI54 method, which is a seven-stage
%method with a fifth-order forward integrator and a fourth-order embedded
%method.
if size(x0,2) > 1
    error("x0 should be pased as column vector!")
end

%Constants for error control
epstol = 0.8;
facmin = 0.1;
facmax = 5;

%Constants for PI error control
% p of not-embedded method = 5, so p+1=6
kI = 0.3/6;
kP = 0.4/6;

h = h0;
t = t0;
x = x0;
X = x;
T = t;

%KuttaNumbers are the inner Xs of the Runge Kutta method
KuttaNumbers = zeros(length(x0),7);
%The butcher tableau is below
%eigth row is best estimate; ninth row is the error measure.
Butcher = [0 1/5 3/10 4/5 8/9 1 1 0 0; ... 
    0 1/5 3/40 44/45 19372/6561 9017/3168 35/384 5179/57600 71/57600; ... 
    0 0 9/40 -56/15 -25360/2187 -355/33 0 0 0; ...
    0 0 0 32/9 64448/6561 46732/5247 500/1113 7571/16695 -71/16695; ...
    0 0 0 0 -212/729 49/176 125/192 393/640 71/1920; ...
    0 0 0 0 0 -5103/18656 -2187/6784 -92097/339200 -17253/339200; ...
    0 0 0 0 0 0 11/84 187/2100 22/525; ...
    0 0 0 0 0 0 0 1/40 -1/40]';

%We need previous r to use the PI controller
rold = 1;
while t < t1
    if t+h >t1
        h = t1-t;
    end
   
    AcceptStep = false;
    KuttaNumbers(:,1) = x;
    % We have an adaptive step, so we go until we find an acceptable step
    while ~AcceptStep
        hButcher = Butcher*h;
        
        
        KuttaTimes = t + hButcher(1:7, 1);
        kuttafs = zeros(length(x0),7); %Kuttafs holds f of each stage so we
        % only have to calculate it once per stage
        kuttafs(:,1) = f(KuttaTimes(1), KuttaNumbers(:,1), params);
        for s = 2:7 % For every stage, calculate the next X 
            KuttaNumbers(:,s) = x;
            for l = 1:s-1
                KuttaNumbers(:,s) = KuttaNumbers(:,s) + hButcher(s, l+1) * kuttafs(:, l);
            end
            kuttafs(:,s) = f(KuttaTimes(s), KuttaNumbers(:,s), params);
        end
        
        e = 0;
        xhat = x;
        for l = 1:7 %Calculate the error by summing over stages
            e = e + hButcher(9, l+1) * kuttafs(:,   l);
            xhat = xhat + hButcher(8, l+1) * kuttafs(:,l);
        end
        r = max(abs(e)./max(abstol, xhat.*reltol));
        
        AcceptStep = (r <= 1.0);
        if AcceptStep % use PI controller
            x = xhat;
            t = t+h;
            X = [X,x];
            T = [T,t];
            h = max(  min((epstol/r)^kI * (rold/r)^kP, facmax), facmin  )*h;
            rold = r;
        else %Use asymptotic controller
            h = max( min( (epstol/r)^(1/6), facmax), facmin) * h;
        end
    end
    
end

X = X';
T = T';
