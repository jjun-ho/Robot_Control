function G = get_Gravity(th)
global Iz1 g L1;

G = Iz1*g/L1*sin(th);
end