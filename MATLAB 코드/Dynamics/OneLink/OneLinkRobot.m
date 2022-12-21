function dydt = OneLinkRobot(t,y)
global I Im m g r Fs Fv tq

ddy = (tq - m*g*r*sin(y(1)) - Fs*sign(y(2)) - Fv*y(2))/(I+Im);

dydt = [y(2); ddy];