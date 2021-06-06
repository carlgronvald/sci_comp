function [x0] = CSTRx0(p)
%CSTRX0 Creates a starting x0 given parameters that lies on the curve
%approximated by the 1D CSTR problem, where all concentrations can be
%approximated by extent of reaction.
T = 273.15;
x0 = [p('CAin') + 1/p('beta')*(p('Tin')-T); p('CBin') + 2/p('beta')*(p('Tin')-T); T]; 
end

