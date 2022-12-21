%% Dynamics 도출
clc
clear all
close all

%% DH 파라미터, Transformation matrix

syms L1 m1 Ic1 Im1 r1 Iz1
syms th1 dth1 
 
I1xx = 0;
I1yy = Iz1;
I1zz = Iz1;
 
d1 = 0; 
a1 = L1; 
al1 = 0; 
 
DH = [th1, d1, a1, al1]; 

T01 = get_HT(th1, d1, a1, al1);

%% 미분 행렬 (rotation)

Qr = [0 -1 0 0 ;
      1  0 0 0 ;
      0  0 0 0 ;
      0  0 0 0];
  
Q1 = Qr;

%% U 행렬

U11 = Q1*T01;       % i=1 j=1

%% Pseudo-inverse

J1(1,1) = 1/2*(-I1xx + I1yy + I1zz);
J1(2,2) = 1/2*( I1xx - I1yy + I1zz);
J1(3,3) = 1/2*( I1xx + I1yy - I1zz);
J1(1,4) = -m1*(L1-r1);
J1(4,1) = J1(1,4);
J1(4,4) = m1;

%% 관성 행렬

D11 = trace(U11*J1*U11.');         % i=1 k=1 j=1 
M11 = simplify(D11);

M = [M11];

%% Coriolis & Centrifugal 행렬

U111 = Q1*Q1*T01;           % i=1 j=1 k=1 

%% n=1 
h111 = trace(U111*J1*U11.');        % i=1 k=1 m=1 j=1   dth1*dth1

h1 = (dth1^2)*h111;

h = simplify([h1]);

%%
n=2;
syms r1 Iz1 
syms g
r11 = [-(L1-r1); 0; 0; 1];
gv = [0 -g 0 0];

G1 = -(m1*gv*U11*r11);          % i=1 j=1
   
G = simplify([G1]);

%%
%[tq1] = M*[ddth1] + h + G
%[ddth1] = inv(M)*([tq1] - h - G)

syms tq1 
DDTH = inv(M)*([tq1] - h - G);

%%
dydt = simplify([dth1; DDTH(1)]);
%%
%matlabFunction(dydt, 'file', 'one_link.m', 'Optimize', false);
