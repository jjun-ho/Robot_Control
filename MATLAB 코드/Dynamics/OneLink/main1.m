%%              
clc
clear all
close all

global m L g u;

dt = 0.005; ft = 5;

q = pi/4; dq = 0;

data = [];

m = 1; L = 1; g = 9.8148; n =1 ; u = 0;

FG = figure('Position', [300 300 600 600], 'Color', [1 1 1])
AX = axes('parent', FG);

hold on
grid on
axis([-1.5 1.5 -1.5 1.5]);

Px = [0 1];
Py = [0 0];

p = plot(Px, Py, '-ob', 'LineWidth', 3)


for cnt=0:dt:ft
    [t,y] = ode45('one_link_ex', [0 dt],[q; dq]);
    
    index = length(y);

    q = y(index,1);
    dq = y(index, 2);
        
    x = L * sin(q);
    y = -L * cos(q);
    Px = [0 x];
    Py = [0 y];
    
    data(n,1) = cnt;
    data(n,2) = q;
    data(n,3) = dq;
    n = n+1;
    
    cmd = sprintf('Time : %2.2f', cnt);
    clc
    disp(cmd)
    
    if rem(n,10) == 0 
        set(p,'XData', Px, 'YData', Py)
   
        drawnow
    end
end
%%
FG2 = figure('Position', [300 300 600 300], 'Color', [1 1 1]);
AX2 = axes('parent', FG2);

plot(data(:,1), data(:,2), 'r')
grid on
hold on
plot(data(:,1), data(:,3), 'b')
%% Regressor
clc
clear all
close all

global  I Im m g r Fs Fv tq;
I = 0.05; Im = 0.05; m = 0.5; g = 9.806; 
r = 0.2; l = 0.5; Fs = 0.1; Fv = 0.1;

dT = 0.002;
q = 0; qdot = 0; 
n = 1;
W1_int = [0, 0, 0, 0];
u = 0;
P = eye(4,4);
theta = [0;0;0;0];

% 10초까지 dT마다 동작 수행
for t = 0 : dT : 5.0
   % 임의의 입력 토크
   tq = sin(1*t) + cos(10*t);    % link1
   
   % ode함수와 상태 공간 방정식을 이용해 q, qdot update
   [st,x] = ode45('OneLinkRobot',[0,dT],[q; qdot]);
   index = size(x); 
   q = x(index(1), 1) ;
   qdot = x(index(1), 2); 
   
   % Regressor 구하기
   W1_int = W1_int + [0, -g*sin(q), -sign(qdot), -qdot]*dT;
   W2 = [qdot, 0, 0, 0];
   Y = W2 - W1_int;
   
   % u값 계산
   u = u + tq*dT;
   
   % Kalman gain
   P = P - P*Y'*inv(1+Y*P*Y')*Y*P;
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
legend('I+Im','mr','Fs','Fv')



