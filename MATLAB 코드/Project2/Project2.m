%% Termproject 02

clear all clc;
%% Regressor

global I1 I2;
global L1 L2;
global Im1 Im2;
global m1 m2;
global g;
global r1 r2;
global Fs1 Fs2;
global Fv1 Fv2;
global tq1 tq2;
 
% 실제 parameter 값
I1 = 0.05; I2 = 0.05; Im1 = 0.05; Im2 = 0.05; m1 = 0.2; m2 = 0.2; g = 9.806; 
Fs1 = 0.1; Fs2 = 0.1; Fv1 = 0.1; Fv2 = 0.1; L1 = 0.5; L2 = 0.5; 
r1 = 0.1; r2 = 0.1; 
 
% 변수 초기화
dT = 0.002;
q1 = 0; q1dot = 0; q2 = 0; q2dot = 0; 
n = 1;
W1_int = zeros(2,10);
u = 0;
P = eye(10,10);
R = eye(2,2);
theta = [0;0;0;0;0;0;0;0;0;0];
 
% 10초까지 dT마다 동작 수행
for t = 0 : dT : 10
   % 임의의 입력 토크
   tq1 = (sin(1*t) + cos(2*t));    % link1
   tq2 = (sin(3*t) + cos(7*t));    % link2
   tq = [tq1 ; tq2];
   
   % ode함수와 상태 공간 방정식을 이용해 q1, q1dot, q2, q2dot update
   [st,x] = ode45('TwoLinkRobot',[0,dT],[q1; q1dot; q2 ; q2dot]);
   index = size(x); 
   q1 = x(index(1), 1) ;
   q1dot = x(index(1), 2); 
   q2 = x(index(1), 3);
   q2dot = x(index(1), 4); 
   
   % Regressor 구하기
   W1_int = W1_int + [ 0 0 0 g*cos(q1) g*cos(q1+q2) 0 -sign(q1dot) 0 -q1dot 0 ; 
                       0 0 -sin(q2)*(q1dot+q2dot)*q1dot 0 g*cos(q1+q2) 0 0 -sign(q2dot) 0 -q2dot] * dT;
   W2 = [ q1dot q1dot+q2dot cos(q2)*(2*q1dot+q2dot) 0 0 0 0 0 0 0 ; 
           0 q1dot+q2dot cos(q2)*q1dot 0 0 q2dot 0 0 0 0 ];
   Y = W2 - W1_int;
   
   % u값 계산
   u = u + tq*dT;
   
   % Kalman gain
   P = P - P*Y'*inv(R+Y*P*Y')*Y*P;
   K = P*Y';
   
   % Theta ipdate
   theta = theta + K*(u - Y*theta);
 
   % Data 저장
   save_time(n,:) = t;
   save_theta(n,:) = theta;
   n=n+1;
 
end
 
% Graph
plot(save_time,save_theta,'DisplayName','save_time')
hold
legend('I1+m2*L2*L2+I1','I2','m2*r2*L1','m1*r1+m2*L2','m2*r2','Im2','Fs1','Fs2','Fv1','Fv2')

