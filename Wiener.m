function [W, tw, dW] = Wiener(t1, intervals, dimensions, realizations, seed)
%WIENER Summary of this function goes here
%   Detailed explanation goes here
if nargin == 4
    rng(seed);
end

dt = T/intervals;
dW = sqrt(dt)*randn(dimensions, intervals, realizations);
W = [zeros(dimensions, 1, realizations) cumsum(dW, 2)];
tw = 0:dt:t1;