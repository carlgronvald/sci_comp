function [params] = CreateParams(varargin)
%CREATEPARAMS Creates a map that can be used to pass parameters to a
%differential equation function
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

