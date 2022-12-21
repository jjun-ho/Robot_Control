%% Termproject 5-1

clc,clear
close all
%% Set Simulation Parameters
    % Draw flag
        flag_Simul      = 1;

        flag_Draw       = 1;
        flag_Draw_Robot = 1;
        flag_Draw_Graph = 1;

    % Global Variable
        global Iz1 Iz2 Iz3 L1 L2 L3 g m1 m2 m3 r1 r2 r3 tq1 tq2 tq3

    % Simulation Parameters
        delta_t     = 0.005;        % [sec]     : Sampling Time - 200Hz
        start_t     = 0.000;        % [sec]     : Start Time
        finish_t    = 5.000;        % [sec]     : End Time 

        g           = 9.8148;       % [m/s^2]   : Gravitational Acceleration
    % Robot Parameters
        m1          = 0.2;          m2          = 0.2;          m3          = 0.2;        % [kg]      : Link Mass
        L1          = 0.5;          L2          = 0.5;          L3          = 0.5;        % [m]       : Link Length
        Iz1         = 0.05;         Iz2         = 0.05;         Iz3         = 0.05;
        r1          = 0.1;          r2          = 0.1;          r3          = 0.1;
        tq1        = 0.000;        tq2        = 0.000;        tq3        = 0.000;      % [Nm]      : Control Torque

        init_q1     = 0.00;         init_q2     = 0.00;         init_q3     = 0.00;       % [rad]     : Init Joint Angle
        init_dq1    = 0.00;         init_dq2    = 0.00;         init_dq3    = 0.00;       % [rad/s]   : Init Angular Velocity
        q1          = init_q1;      q2          = init_q2;      q3          = init_q3;    % [rad]     : Joint Angle
        dq1         = init_dq1;     dq2         = init_dq2;     dq3         = init_dq3;   % [rad/s]   : Angular Velocity
        final_q1    = 90*pi/180;    final_q2    = 90*pi/180;    final_q3    = 90*pi/180;  % [rad]     : Final Joint Angle
        final_dq1   = 30*pi/180;    final_dq2   = 30*pi/180;    final_dq3   = 30*pi/180;  % [rad/s]   : Final Angular Velocity

        init_q      = [  init_q1;    init_q2;    init_q3];
        init_dq     = [ init_dq1;   init_dq2;   init_dq3];
        q           = [       q1;         q2;         q3];
        dq          = [      dq1;        dq2;        dq3];
        final_q     = [ final_q1;   final_q2;   final_q3];
        final_dq    = [final_dq1;  final_dq2;  final_dq3];

    % Target Pos/Vel/Acc Parameters
        q1_d        = init_q1;     q2_d         = init_q2;     q3_d         = init_q3;   % [rad]     : Target Joint Angle
        dq1_d       = init_dq1;    dq2_d        = init_dq2;    dq3_d        = init_dq3;  % [rad/s]   : Target Angular Velocity
        ddq1_d      = 0;           ddq2_d       = 0;           ddq3_d       = 0;         % [rad/s^2] : Target Angular Acceleration

        q_d         = [       q1_d;         q2_d;         q3_d];
        dq_d        = [      dq1_d;        dq2_d;        dq3_d];
        ddq_d       = [     ddq1_d;       ddq2_d;       ddq3_d];

    % Controller Gain
        Wn          = 20;                                         % [rad/s]   : natural frequency
        Kp          = [Wn^2;        Wn^2;       Wn^2        ];    % [Nm/rad]  : Propotional Gain
        zeta        = 1;                                          % Damping ratio(감쇠비)
        Kv          = [2*zeta*Wn;   2*zeta*Wn;  2*zeta*Wn   ];    % [Nm*s/rad]: Derivative Gain
        Ki          = [400;         400;        400         ];
        error       = zeros(3,1);
        errorSum    = zeros(3,1);
        u           = zeros(3,1);
%% Simulation
    if(flag_Simul == 1)
        % Simulation
            n = 1;
            errorSum = 0;
            for (time = start_t:delta_t:finish_t)
                % Set Target Trajectory
                for (i = 1:3)
                    if(time<1)
                        q_d(i)     = init_q(i);
                        dq_d(i)    = init_dq(i);
                        ddq_d(i)   = 0.0;
                    else
                        if(q_d(i) < final_q(i))
                            q_d(i) = q_d(i) + (final_dq(i))*delta_t;
                        elseif(q_d(i) > final_q(i) + 0.0001)
                            q_d(i) = q_d(i) - (final_dq(i))*delta_t;
                        else 
                            q_d(i) = final_q(i);
                        end
                        dq_d(i)    = (q_d(i)  - simul_q_d(i,n-1)) / delta_t;
                        ddq_d(i)   = (dq_d(i) - simul_dq_d(i,n-1))/ delta_t;
                    end

                    % Get Dynamics
                    G           = get_Gravity3(q(1),q(2),q(3));
                    % Controller
                    error       = q_d - q;
                    errorSum    = errorSum + error*delta_t;
                    u           = ddq_d + Kv(i)*(dq_d - dq) + Kp(i)*(q_d - q) + Ki(i)*errorSum;
                    I = get_Inertia3(q(1),q(2),q(3));
                    tq_ctrl    = I*u + G * 0.8;

                end
                

                % Robot Model  
                    % Inverse Dynamics
                        tq1      = tq_ctrl(1);
                        tq2      = tq_ctrl(2);
                        tq3      = tq_ctrl(3);
                        [t,y]   = ode45('three_link',[0 delta_t],[q(1); dq(1); q(2); dq(2); q(3); dq(3)]);
                        index   = length(y);
                        q(1)    = y(index,1);
                        dq(1)   = y(index,2);
                        q(2)    = y(index,3);
                        dq(2)   = y(index,4);
                        q(3)    = y(index,5);
                        dq(3)   = y(index,6);

                    % Save Data
                        simul_time(n)   = time;     % [sec]
                        simul_q(:,n)    = q;        % [rad]
                        simul_dq(:,n)   = dq;       % [rad/s]
                        simul_q_d(:,n)  = q_d;      % [rad]
                        simul_dq_d(:,n) = dq_d;     % [rad/s]
                        n               = n + 1;

                    % Process
                        cmd = sprintf("Loading(%2.0f%%)", time/finish_t*100);
                        clc
                        disp(cmd)
            end
    end
%% Simulation Result Graph
    if(flag_Draw == 1)
        font_size_label     = 20;
        font_size_title     = 25;
        linewidth_current   = 3;
        linewidth_target    = 5;

        if(flag_Draw_Robot == 1)
            % Draw Robot
                init_x1      = L1*cos(init_q(1));
                init_y1      = L1*sin(init_q(1));

                init_x2      = L2*cos(init_q(1)+init_q(2));
                init_y2      = L2*sin(init_q(1)+init_q(2));

                init_x3      = L3*cos(init_q(1)+init_q(2)+init_q(3));
                init_y3      = L3*sin(init_q(1)+init_q(2)+init_q(3));
            
                FG1 = figure('Position',[200 300 500 450],'Color',[1 1 1]);
                    Ax = axes('Parent',FG1); hold on

                    p1 = plot([0 0]              ,[init_x1 init_y1],'-ob','LineWidth',linewidth_current);
                    p2 = plot([init_x1 init_y1]  ,[init_x2 init_y2],'-or','LineWidth',linewidth_current);
                    p3 = plot([init_x2 init_y2]  ,[init_x3 init_y3],'-og','LineWidth',linewidth_current);

                    axis([-1.5 1.5 -1.5 1.5]);
                    grid on
                xlabel('X-axis (m)',    'FontSize',font_size_label)
                ylabel('Y-axis (m)',    'FontSize',font_size_label)
                title('3-DOF Robot',    'FontSize',font_size_title)

                n = 1;
                for (time = start_t:delta_t:finish_t)
                    q(1) = simul_q(1,n);
                    x1 = L1*cos(q(1));
                    y1 = L1*sin(q(1));
                    Px1 = [0 x1];
                    Py1 = [0 y1];
                    
                    q(2) = simul_q(2,n);
                    x2 = L2*cos(q(1)+q(2));
                    y2 = L2*sin(q(1)+q(2));
                    Px2 = [x1 x1+x2];
                    Py2 = [y1 y1+y2];
                    
                    q(3) = simul_q(3,n);
                    x3 = L3*cos(q(1)+q(2)+q(3));
                    y3 = L3*sin(q(1)+q(2)+q(3));
                    Px3 = [x1+x2 x1+x2+x3];
                    Py3 = [y1+y2 y1+y2+y3];

                    set(p1, 'XData',Px1,'YData',Py1)
                    set(p2, 'XData',Px2,'YData',Py2)
                    set(p3, 'XData',Px3,'YData',Py3)
                    drawnow
                    n = n +1;
                    
                    cmd = sprintf("Time : %2.2f",time);
                    clc
                    disp(cmd)
                end
        end
        if(flag_Draw_Graph == 1)
            % Draw Angle
                FG2 = figure('Position',[800 700 600 300],'Color',[1 1 1]);
                    plot(simul_time,simul_q_d * 180/pi, ':k','LineWidth',linewidth_target);
                    hold on;
                    plot(simul_time,simul_q * 180/pi, ' r','LineWidth',linewidth_current);
                    hold on;

                    legend('Desired','Current')
                    axis([start_t finish_t 0 120]);
                    xticks([start_t:1:finish_t])
                    yticks([0:45:90])
                    grid on
                xlabel('time (s)',                      'FontSize',font_size_label)
                ylabel('Angle (deg)',                   'FontSize',font_size_label)
                title('Joint Space PID CTM Controller', 'FontSize',font_size_title)

            % Draw Angular Velocity
                FG3 = figure('Position',[800 100 600 300],'Color',[1 1 1]);
                    plot(simul_time,simul_dq_d * 180/pi, ':k','LineWidth',linewidth_target);
                    hold on;
                    plot(simul_time,simul_dq * 180/pi, ' r','LineWidth',linewidth_current);
                    hold on;

                    legend('Desired','Current')
                    axis([start_t finish_t -90 90]);
                    xticks([start_t:1:finish_t])
                    yticks([-60 0 30 60])
                    grid on
                xlabel('time (s)',                      'FontSize',font_size_label)
                ylabel('Angular Velocity (deg/s)',      'FontSize',font_size_label)
                title('Joint Space PID CTM Controller', 'FontSize',font_size_title)
        end
    end