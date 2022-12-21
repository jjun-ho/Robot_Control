function Gout = get_Gravity3(th1, th2, th3)
global g L1 L2 L3 r1 r2 r3 m1 m2 m3;

G1 = g*m2*(L1*cos(th1) + L2*cos(th1)*cos(th2) - L2*sin(th1)*sin(th2)) - g*m3*(sin(th1)*(L2*sin(th2) + L3*cos(th2)*sin(th3) + L3*cos(th3)*sin(th2)) - cos(th1)*(L2*cos(th2) + L3*cos(th2)*cos(th3) - L3*sin(th2)*sin(th3))) - g*m2*(L2 - r2)*(cos(th1)*cos(th2) - sin(th1)*sin(th2)) + L1*g*m1*cos(th1) - g*m1*cos(th1)*(L1 - r1) - g*m3*(cos(th1)*(cos(th2)*cos(th3) - sin(th2)*sin(th3)) - sin(th1)*(cos(th2)*sin(th3) + cos(th3)*sin(th2)))*(L3 - r3);
G2 = g*m2*(L2*cos(th1)*cos(th2) - L2*sin(th1)*sin(th2)) - g*m3*(sin(th1)*(L2*sin(th2) + L3*cos(th2)*sin(th3) + L3*cos(th3)*sin(th2)) - cos(th1)*(L2*cos(th2) + L3*cos(th2)*cos(th3) - L3*sin(th2)*sin(th3))) - g*m2*(L2 - r2)*(cos(th1)*cos(th2) - sin(th1)*sin(th2)) - g*m3*(cos(th1)*(cos(th2)*cos(th3) - sin(th2)*sin(th3)) - sin(th1)*(cos(th2)*sin(th3) + cos(th3)*sin(th2)))*(L3 - r3);
G3 = - g*m3*(sin(th1)*(L2*sin(th2) + L3*cos(th2)*sin(th3) + L3*cos(th3)*sin(th2)) - cos(th1)*(L2*cos(th2) + L3*cos(th2)*cos(th3) - L3*sin(th2)*sin(th3))) - g*m3*(cos(th1)*(cos(th2)*cos(th3) - sin(th2)*sin(th3)) - sin(th1)*(cos(th2)*sin(th3) + cos(th3)*sin(th2)))*(L3 - r3);

Gout = [G1; G2; G3];
end