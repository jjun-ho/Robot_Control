function dxdt = TwoLinkRobot(t,x)

global I1 I2;
global L1 L2;
global Im1 Im2;
global m1 m2;
global g;
global r1 r2;
global Fs1 Fs2;
global Fv1 Fv2;
global tq1 tq2;

q1 = x(1);
q2 = x(3);
v1 = x(2);
v2 = x(4);

M_11 = I1 + I2 + m2*L1*L1 + 2*m2*r2*L1*cos(q2) + Im1;
M_12 = I2 + m2*r2*L1*cos(q2);
M_21 = I2 + m2*r2*L1*cos(q2);
M_22 = I2 + Im2;
C_11 = -m2*r2*L1*sin(q2)*v2;
C_12 = -m2*r2*L1*sin(q2)*(v1+v2);
C_21 = m2*r2*L1*sin(q2)*v1;
C_22 = 0;
g_1 = -m1*r1*g*cos(q1) - m2*L1*g*cos(q1) - m2*r2*g*cos(q1+q2);
g_2 = -m2*r2*g*cos(q1+q2);
d_1 = Fs1*sign(v1) + Fv1*v1;
d_2 = Fs2*sign(v2) + Fv2*v2;
 
C = [C_11 C_12 ; C_21 C_22];
G = [g_1 ; g_2];
D = [d_1 ; d_2];
R = C*[v1 ; v2] + G + D;
 
tq = [tq1 ; tq2] - R;
M = [M_11, M_12 ; M_21, M_22];
vdot = M\tq;
dxdt = [v1 ; vdot(1) ; v2 ; vdot(2)];
end
