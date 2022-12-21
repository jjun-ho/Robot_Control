function [X] = get_Kinematics3(th1, th2, th3)
    global L1 L2 L3

    x = L1*cos(th1) + L2*cos(th1+th2) + L3*cos(th1+th2+th3);
    y = L1*sin(th1) + L2*sin(th1+th2) + L3*sin(th1+th2+th3);

    X = [x;y];
end