function [W, tw, dW] = Wiener(t0, t1, intervals, dimensions, realizations, seed)
%WIENER realizes a standard Wiener process
% t0 - starting time
% t1 - ending time
% intervals - number of steps between t0 and t1
% dimensions - How many dimensions do we need
% realizations - How many times do we realize it
% seed - seed for the random number generator
if nargin == 6
    rng(seed);
end

dt = (t1-t0)/intervals;
dW = sqrt(dt)*randn(dimensions, intervals, realizations);
W = [zeros(dimensions, 1, realizations) cumsum(dW, 2)];
tw = 0:dt:t1;