function [params] = CreateParams(varargin)
%CREATEPARAMS Summary of this function goes here
%   Detailed explanation goes here
keys = [];
values = [];
for i=1:nargin
    if mod(i, 2) == 1
        keys = [keys; varargin(i)];
    else
        values = [values; varargin(i)];
    end
end
params = containers.Map(keys, values);

