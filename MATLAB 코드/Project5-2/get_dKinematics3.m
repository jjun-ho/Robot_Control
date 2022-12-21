function [dXout] = get_dKinematics3(th1, th2, th3)
    global L1 L2 L3

    x = -L1*sin(th1)*(th1) - L2*sin(th1+th2)*(th1+th2) - L3*sin(th1+th2+th3)*(th1+th2+th3);
    y = L1*cos(th1)*(th1) + L2*cos(th1+th2)*(th1+th2) + L3*cos(th1+th2+th3)*(th1+th2+th3);

    dXout = [x;y];
end