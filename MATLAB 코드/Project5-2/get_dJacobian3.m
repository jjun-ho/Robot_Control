function [dJout] = get_dJacobian3(th1, th2, th3)
global L1 L2 L3

J11 = -L2*cos(th1+th2)*(th1+th2)-L1*cos(th1)*(th1)-L3*cos(th1+th2+th3)*(th1+th2+th3);
J12 = -L2*cos(th1+th2)*(th1+th2)-L3*cos(th1+th2+th3)*(th1+th2+th3);
J13 = -L3*cos(th1+th2+th3)*(th1+th2+th3);
J21 = -L2*sin(th1+th2)*(th1+th2)-L1*sin(th1)*(th1)-L3*sin(th1+th2+th3)*(th1+th2+th3);
J22 = -L2*sin(th1+th2)*(th1+th2)-L3*sin(th1+th2+th3)*(th1+th2+th3);
J23 = -L3*sin(th1+th2+th3)*(th1+th2+th3);

dJout = [J11, J12, J13;
        J21, J22,  J23];
end