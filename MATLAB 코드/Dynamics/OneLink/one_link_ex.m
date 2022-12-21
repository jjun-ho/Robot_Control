function [dydt] = one_link_ex(t, y)
global m L g u;

dydt = [y(2); -g/L*sin(y(1)) + u/m/L/L];
