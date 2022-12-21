%% Dynamics 도출
clc
clear all
close all

%% DH 파라미터, Transformation matrix
syms L1 L2 L3 m1 m2 m3 Ic1 Ic2 Ic3 Im1 Im2 Im3 r1 r2 r3 Iz1 Iz2 Iz3
syms th1 th2 th3 dth1 dth2 dth3 ddth1 ddth2 ddth3

I1xx = 0;
I1yy = Iz1;
I1zz = Iz1;

I2xx = 0;
I2yy = Iz2;
I2zz = Iz2;

I3xx = 0;
I3yy = Iz3;
I3zz = Iz3;
 
d1 = 0; d2 = 0; d3 = 0;
a1 = L1; a2 = L2; a3 = L3;
al1 = 0; al2 = 0; al3 = 0;
 
DH = [th1, d1, a1, al1 ; 
      th2, d2, a2, al2 ;
      th3, d3, a3, al3];

T01 = get_HT(th1, d1, a1, al1);
T12 = get_HT(th2, d2, a2, al2);
T23 = get_HT(th3, d3, a3, al3);
T02 = T01 * T12; 
T13 = T12 * T23;
T03 = T01 * T12 * T23;

%% 미분 행렬 (rotation)
Qr = [0 -1 0 0 ;
      1  0 0 0 ;
      0  0 0 0 ;
      0  0 0 0];
  
Q1 = Qr;
Q2 = Qr;
Q3 = Qr;

%% U 행렬

U11 = Q1*T01;       % i=1 j=1
U12 = zeros(4,4);   % i=1 j=2
U13 = zeros(4,4);   % i=1 j=3

U21 = Q1*T02;       % i=2 j=1
U22 = T01*Q2*T12;   % i=2 j=2
U23 = zeros(4,4);   % i=2 j=3

U31 = Q3*T03;       % i=3 j=1
U32 = T01*Q2*T13;   % i=3 j=2
U33 = T02*Q3*T23;   % i=3 j=3

%% Pseudo-inverse

J1(1,1) = 1/2*(-I1xx + I1yy + I1zz);
J1(2,2) = 1/2*( I1xx - I1yy + I1zz);
J1(3,3) = 1/2*( I1xx + I1yy - I1zz);
J1(1,4) = -m1*(L1-r1);
J1(4,1) = J1(1,4);
J1(4,4) = m1;

J2(1,1) = 1/2*(-I2xx + I2yy + I2zz);
J2(2,2) = 1/2*( I2xx - I2yy + I2zz);
J2(3,3) = 1/2*( I2xx + I2yy - I2zz);
J2(1,4) = -m2*(L2-r2);
J2(4,1) = J2(1,4);
J2(4,4) = m2;

J3(1,1) = 1/2*(-I3xx + I3yy + I3zz);
J3(2,2) = 1/2*( I3xx - I3yy + I3zz);
J3(3,3) = 1/2*( I3xx + I3yy - I3zz);
J3(1,4) = -m3*(L3-r3);
J3(4,1) = J3(1,4);
J3(4,4) = m3;

%% 관성 행렬

D11 = trace(U11*J1*U11.')...       % i=1 k=1 j=1 
    + trace(U21*J2*U21.')...       % i=1 k=1 j=2
    + trace(U31*J3*U31.');         % i=1 k=1 j=3
M11 = simplify(D11);

D12 = trace(U22*J2*U21.')...       % i=1 k=2 j=2 
    + trace(U32*J3*U31.');         % i=1 k=2 j=3
M12 = simplify(D12);

D13 = trace(U33*J3*U31.');         % i=1 k=3 j=3 
M13 = simplify(D13);

M21 = M12;                         % D21 = trace(U21*J2*U22.')

D22 = trace(U22*J2*U22.')...       % i=2 k=2 j=2 
    + trace(U32*J3*U32.');         % i=2 k=2 j=3
M22 = simplify(D22);

D23 = trace(U33*J3*U32.');         % i=2 k=3 j=3 
M23 = simplify(D23);

M31 = M13;

M32 = M23;

D33 = trace(U33*J3*U33.');          % i=3 k=3 j=3 
M33 = simplify(D33);

M = [M11 M12 M13; M21 M22 M23; M31 M32 M33];

%% Coriolis & Centrifugal 행렬

n = 3;

for i=1:n
    for j=1:n
        for k=1:n
            cmd = sprintf("U%d%d%d = get_dUdq3(i,j,k,T01,T12,T23,T02,T13,T03,Q1,Q2,Q3);",i,j,k);
            eval(cmd);
        end
    end
end

%% n=3
h111 = trace(U111*J1*U11.')...      % i=1 k=1 m=1 j=1
     + trace(U211*J2*U21.')...      % i=1 k=1 m=1 j=2   
     + trace(U311*J3*U31.');        % i=1 k=1 m=1 j=3   dth1*dth1
h112 = trace(U212*J2*U21.')...      % i=1 k=1 m=2 j=2   
     + trace(U212*J2*U21.');        % i=1 k=1 m=2 j=3   dth1*dth2
h113 = trace(U312*J3*U31.');        % i=1 k=1 m=3 j=3   dth1*dth3

h121 = trace(U221*J2*U21.')...      % i=1 k=2 m=1 j=2   
     + trace(U321*J3*U31.');        % i=1 k=2 m=1 j=3   dth2*dth1
h122 = trace(U222*J2*U21.')...      % i=1 k=2 m=2 j=2   
     + trace(U322*J3*U31.');        % i=1 k=2 m=2 j=3   dth2*dth2
h123 = trace(U323*J3*U31.');        % i=1 k=2 m=3 j=3   dth2*dth3

h131 = trace(U331*J3*U31.');        % i=1 k=3 m=1 j=3   dth3*dth1
h132 = trace(U332*J3*U31.');        % i=1 k=3 m=2 j=3   dth3*dth2
h133 = trace(U333*J3*U31.');        % i=1 k=3 m=3 j=3   dth3*dth3

h1 = (dth1^2)*h111 + (dth1*dth2)*(h112+h121) + (dth1*dth3)*(h113+h131) + (dth2^2)*h122 + (dth2*dth3)*(h123+h132) + (dth3^2)*h133;

h211 = trace(U111*J1*U12.')...      % i=2 k=1 m=1 j=1
     + trace(U211*J2*U22.')...      % i=2 k=1 m=1 j=2   
     + trace(U311*J3*U32.');        % i=2 k=1 m=1 j=3   dth1*dth1
h212 = trace(U212*J2*U22.')...      % i=2 k=1 m=2 j=2   
     + trace(U212*J2*U22.');        % i=2 k=1 m=2 j=3   dth1*dth2
h213 = trace(U312*J3*U32.');        % i=2 k=1 m=3 j=3   dth1*dth3

h221 = trace(U221*J2*U22.')...      % i=2 k=2 m=1 j=2   
     + trace(U321*J3*U32.');        % i=2 k=2 m=1 j=3   dth2*dth1
h222 = trace(U222*J2*U22.')...      % i=2 k=2 m=2 j=2   
     + trace(U322*J3*U32.');        % i=2 k=2 m=2 j=3   dth2*dth2
h223 = trace(U323*J3*U32.');        % i=2 k=2 m=3 j=3   dth2*dth3

h231 = trace(U331*J3*U32.');        % i=2 k=3 m=1 j=3   dth3*dth1
h232 = trace(U332*J3*U32.');        % i=2 k=3 m=2 j=3   dth3*dth2
h233 = trace(U333*J3*U32.');        % i=2 k=3 m=3 j=3   dth3*dth3

h2 = (dth1^2)*h211 + (dth1*dth2)*(h212+h221) + (dth1*dth3)*(h213+h231) + (dth2^2)*h222 + (dth2*dth3)*(h223+h232) + (dth3^2)*h233;

h311 = trace(U111*J1*U13.')...      % i=3 k=1 m=1 j=1
     + trace(U211*J2*U23.')...      % i=3 k=1 m=1 j=2   
     + trace(U311*J3*U33.');        % i=3 k=1 m=1 j=3   dth1*dth1
h312 = trace(U212*J2*U23.')...      % i=3 k=1 m=2 j=2   
     + trace(U212*J2*U23.');        % i=3 k=1 m=2 j=3   dth1*dth2
h313 = trace(U312*J3*U33.');        % i=3 k=1 m=3 j=3   dth1*dth3

h321 = trace(U221*J2*U23.')...      % i=3 k=2 m=1 j=2   
     + trace(U321*J3*U33.');        % i=3 k=2 m=1 j=3   dth2*dth1
h322 = trace(U222*J2*U23.')...      % i=3 k=2 m=2 j=2   
     + trace(U322*J3*U33.');        % i=3 k=2 m=2 j=3   dth2*dth2
h323 = trace(U323*J3*U33.');        % i=3 k=2 m=3 j=3   dth2*dth3

h331 = trace(U331*J3*U33.');        % i=3 k=3 m=1 j=3   dth3*dth1
h332 = trace(U332*J3*U33.');        % i=3 k=3 m=2 j=3   dth3*dth2
h333 = trace(U333*J3*U33.');        % i=3 k=3 m=3 j=3   dth3*dth3

h3 = (dth1^2)*h311 + (dth1*dth2)*(h312+h321) + (dth1*dth3)*(h313+h331) + (dth2^2)*h322 + (dth2*dth3)*(h323+h332) + (dth3^2)*h333;

h = simplify([h1; h2; h3]);

%%
n=3;
syms r1 r2 r3 Iz1 Iz2 Iz3
syms g
r11 = [-(L1-r1); 0; 0; 1];
r22 = [-(L2-r2); 0; 0; 1];
r33 = [-(L3-r3); 0; 0; 1];
gv = [0 -g 0 0];

G1 = -(m1*gv*U11*r11 +...       % i=1 j=1
       m2*gv*U21*r22 +...       % i=1 j=2
       m3*gv*U32*r33);          % i=1 j=3
   
G2 = -(m2*gv*U22*r22 +...       % i=2 j=2
       m3*gv*U32*r33);          % i=2 j=3 

G3 = -(m3*gv*U32*r33);          % i=3 j=3

G = simplify([G1; G2; G3]);

%%

%[tq1; tq2; tq3] = M*[ddth1;ddth2;ddth3] + h + G
%[ddth1;ddth2; tq3] = inv(M)*([tq1; tq2; tq3] - h - G)

syms tq1 tq2 tq3
DDTH = inv(M)*([tq1; tq2; tq3] - h - G);
%%
dydt = simplify([dth1; DDTH(1); dth2; DDTH(2); dth3; DDTH(3)]);
%%
%matlabFunction(dydt, 'file', 'three_link.m', 'Optimize', false);
