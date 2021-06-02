function [x0] = CSTRx0(p)
%CSTRX0 Summary of this function goes here
%   Detailed explanation goes here
T = 273.15;
x0 = [p('CAin') + 1/p('beta')*(p('Tin')-T); p('CBin') + 2/p('beta')*(p('Tin')-T); T]; 
end

