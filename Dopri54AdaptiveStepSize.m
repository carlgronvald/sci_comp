function [X,T] = Dopri54AdaptiveStepSize(x0, f, h0, t0, t1, abstol, reltol, params)
%DOPRI54 Summary of this function goes here
%   Detailed explanation goes here
if size(x0,2) > 1
    error("x0 should be pased as column vector!")
end

%Constants for asymptotic error control
epstol = 0.8;
facmin = 0.1;
facmax = 5;

%Constants for PI error control
% p of not-embedded method = 4, so p+1=5
kI = 0.3/5;
kP = 0.4/5;

h = h0;
t = t0;
x = x0;
X = x;
T = t;

KuttaConstants = zeros(length(x0),7);
KuttaTimes = zeros(1,7);
%eigth row is best estimate; ninth row is error estimate
Butcher = [0 1/5 3/10 4/5 8/9 1 1 0 0; ... 
    0 1/5 3/40 44/45 19372/6561 9017/3168 35/384 5179/57600 71/57600; ... 
    0 0 9/40 -56/15 -25360/2187 -355/33 0 0 0; ...
    0 0 0 32/9 64448/6561 46732/5247 500/1113 7571/16695 -71/16695; ...
    0 0 0 0 -212/729 49/176 125/192 393/640 71/1920; ...
    0 0 0 0 0 -5103/18656 -2187/6784 -92097/339200 -17253/339200; ...
    0 0 0 0 0 0 11/84 187/2100 22/525; ...
    0 0 0 0 0 0 0 1/40 -1/40]';

rold = 1;
while t < t1
    if t+h >t1
        h = t1-t;
    end
   
    AcceptStep = false;
    KuttaNumbers(:,1) = x;
    while ~AcceptStep

        %disp("Not accepting step");
        hButcher = Butcher*h;
        
        
       % KuttaConstants(:,1) = X(:,i);
        %KuttaTimes = T(i) + hButcher(1:4, 1);
        %for s = 2:4
        %    KuttaConstants(:,s) = X(:,i) + hButcher(s, 2:s) * f(KuttaTimes(1:s-1),KuttaConstants(1,1:s-1),params)';
        %end

        %X(:,i+1) = X(:,i) +  hButcher(5,2:5) * f(KuttaTimes, KuttaConstants, params)';
        %T(i+1) = T(i)+h;
        
        KuttaTimes = t + hButcher(1:7, 1);
        for s = 2:7
            KuttaNumbers(:,s) = x + hButcher(s, 2:s) * f(KuttaTimes(1:s-1),KuttaNumbers(1,1:s-1),params)';
        end
        
    	e = hButcher(9,2:8) * f(KuttaTimes, KuttaNumbers, params)';
        xhat = x + hButcher(8,2:8) * f(KuttaTimes, KuttaNumbers, params)';
        r = max(abs(e)./max(abstol, xhat.*reltol));
        
        AcceptStep = (r <= 1.0);
        if AcceptStep % use PI controller
            x = xhat;
            t = t+h;
            X = [X;x];
            T = [T;t];
            h = max(  min((epstol/r)^kI * (rold/r)^kP, facmax), facmin  )*h;
            rold = r;
        else
            h = max( min( (epstol/r)^(1/5), facmax), facmin) * h;
        end
    end
    
end

