function [dJout] = get_dJacobian(th1,th2)
    global L1 L2
     
    J11 = -L1*cos(th1)*(th1)-L2*cos(th1+th2)*(th1+th2);    
    J12 = -L2*cos(th1+th2)*(th1+th2);
    J21 = -L1*sin(th1)*(th1)-L2*sin(th1+th2)*(th1+th2);
    J22 = -L2*sin(th1+th2)*(th1+th2);
    
    dJout = [J11, J12; J21, J22];
end
    