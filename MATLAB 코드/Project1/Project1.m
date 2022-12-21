%% Termproject 01

clear all
close all
global Iz1 Iz2 Iz3 L1 L2 L3 g m1 m2 m3 r1 r2 r3 tq1 tq2 tq3

L1 = 0.5;   L2 = 0.5;   L3 = 0.5;
r1 = 0.1;   r2 = 0.1;   r3 = 0.1;
m1 = 0.2;   m2 = 0.2;   m3 = 0.2;
Iz1=0.05;   Iz2=0.05;   Iz3=0.05;

g = 9.806;

dt = 0.005; ft = 5;

%초기각도
q1 = -pi/2; dq1 = 0;
q2 =  pi/4; dq2 = 0;
q3 =  pi/4; dq3 = 0;

data = [];

n = 1;

FG = figure('Position',[200 300 500 450], 'Color',[1 1 1])
AX = axes('Parent',FG);

hold on
grid on
axis([-1.5 1.5 -1.5 1.5]);

x1 = L1*cos(q1);
y1 = L1*sin(q1);
Px1 = [0 x1];
Py1 = [0 y1];

x2 = L2*cos(q1+q2);
y2 = L2*sin(q1+q2);
Px2 = [x1 x1+x2];
Py2 = [y1 y1+y2];

x3 = L3*cos(q1+q2+q3);
y3 = L3*sin(q1+q2+q3);
Px3 = [x1+x2 x1+x2+x3];
Py3 = [y1+y2 y1+y2+y3];

p1 = plot(Px1, Py1,'-ob','Linewidth',3)
p2 = plot(Px2, Py2,'-or','Linewidth',3)
p3 = plot(Px3, Py3,'-og','Linewidth',3)


for cnt=0:dt:ft
 
    tq1 = 0.0;
    tq2 = 0.0;
    tq3 = 0.0;

    [t,y] = ode45('three_link',[0 dt],[q1; dq1; q2; dq2; q3; dq3]);

    index = length(y);

    q1  = y(index,1);
    dq1 = y(index,2);
    q2  = y(index,3);
    dq2 = y(index,4);
    q3  = y(index,5);
    dq3 = y(index,6);

    x1 = L1*cos(q1);
    y1 = L1*sin(q1);
    Px1 = [0 x1];
    Py1 = [0 y1];
    
    x2 = L2*cos(q1+q2);
    y2 = L2*sin(q1+q2);
    Px2 = [x1 x1+x2];
    Py2 = [y1 y1+y2];
    
    x3 = L3*cos(q1+q2+q3);
    y3 = L3*sin(q1+q2+q3);
    Px3 = [x1+x2 x1+x2+x3];
    Py3 = [y1+y2 y1+y2+y3];

    data(n,1) = cnt;
    data(n,2) = q3;
    data(n,3) = dq3;
    
    n = n + 1;

    cmd = sprintf("Time : %2.2f",cnt);
    clc
    disp(cmd)

    if rem(n,10) == 0
        set(p1, 'XData',Px1,'YData',Py1)
        set(p2, 'XData',Px2,'YData',Py2)
        set(p3, 'XData',Px3,'YData',Py3)
        drawnow
    end
end
%%
FG2 = figure('Position', [800 400 600 300], 'Color', [1 1 1])
AX2 = axes('parent', FG2);

plot(data(:,1), data(:,2), 'r')
grid on
hold on
plot(data(:,1), data(:,3), 'b')

legend('q(rad)', 'dq(rad/s)');

xlabel('time (s)', 'fontsize', 20)
ylabel('Y-axis', 'fontsize', 20) 
title('3-DOF Robot and free-fall Simulation', 'fontsize', 25)

