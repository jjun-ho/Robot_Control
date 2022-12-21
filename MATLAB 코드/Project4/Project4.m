%% Termproject 04

clc, clear all
close all
%% Set Simulation Parameters
    % Draw Flag
        flag_Simul = 1;
        flag_Draw = 1;
        flag_Draw_Robot = 1;
        flag_Draw_Graph = 1;
    % Global Varial
        global Iz1 Iz2 L1 L2 g m1 m2 r1 r2 tq1 tq2;
    % Simulation Parameters
        delta_t     = 0.005;        % [sec]     : Sampling Time
        start_t     = 0.000;        % [sec]     : Start Time
        finish_t    = 5.000;        % [sec]     : End Time
        
        g           = 9.8148;       % [m/s^2]   : Gravitational Acceleration
    % Robot Parameters
        m1  = 0.20; m2 = 0.20;      % [kg]      : Link mass
        L1  = 0.50; L2 = 0.50;      % [m]       : Link Length
        r1  = 0.10; r2 = 0.10;      % [m]       : Center of Length
        Iz1 = 0.10; Iz2 = 0.05;     % [kgm^2]   : Link Inertia
                
        init_q1  = -pi/2;    init_q2 = pi/2;     % [rad]     : Init Joint Angle
        init_dq1 = 0.00;     init_dq2= 0.00;     % [rad/s]   : Init Angular Velocity
        q           = [init_q1
                       init_q2];                 % [rad]     : Init Joint Angle
        dq          = [init_dq1
                       init_dq2];                % [rad/s]   : Init Angular Velocity

        init_X  = get_Kinematics(q(1),q(2));     % [m]       : Init End-Effector Position
        X       = init_X;                        % [m]       : Current Position
        dX      = [0;0];                         % [m]       : Current Velocity
        X_d     = init_X;                        % [m]       : Current End-Effector Position  
        dX_d    = [0;0];                         % [m/s]     : Current End-Effector Velocity      
        ddX_d   = [0;0];                         % [m/s^2]   : Current End-Effector Acceleration
                  
        tq1 = 0.0;                               % [Nm]      : Control Torque 1
        tq2 = 0.0;                               % [Nm]      : Control Torque 2
        tq  = [ tq1
                tq2];                            % [Nm]      : Control Torque   
        % Controller Parameters
            Wn          = 50;                          % [rad/s]   : natural frequency
            Kp          = [Wn^2;        Wn^2];         % [Nm/rad]  : Propotional Gain
            zeta        = 1;                           % Damping ratio(감쇠비)
            Kv          = [2*zeta*Wn;   2*zeta*Wn];    % [Nm*s/rad]: Derivative Gain
            Ki          = [800;         800;];
            error       = zeros(2,1);
            errorSum    = zeros(2,1);
            u           = zeros(2,1);
%% Simulation
    if(flag_Simul == 1)
        %Simulation
            n = 1;
            errorSum =0;
            sin_t = 0;
            for (time = start_t:delta_t:finish_t)
                % Set Target Trajectory
                for(i = 1:2)
                    if(time<1)
                        X_d         = init_X;
                        dX_d        = [0;0];
                        ddX_d       = [0;0];
                    elseif(time < 2.0)
                        X_d(1)      = init_X(1);
                        if(X_d(2) < init_X(2) + 0.1)    % TargetPos
                            X_d(2) = X_d(2) + (0.1/0.5)*delta_t; % TargetVelocity
                        else
                            X_d(2) = init_X(2) + 0.1;
                        end
                        dX_d  = (X_d  - [simul_X_d_x(n-1); simul_X_d_y(n-1)])./delta_t;
                        ddX_d = (dX_d  - [simul_dX_d_x(n-1); simul_dX_d_y(n-1)])./delta_t;
                    else
                        X_d         = [ 0.1 * sin((2*pi*sin_t)/2) + init_X(1);
                                        0.1 * cos((2*pi*sin_t)/2) + init_X(2)];
                        sin_t = sin_t + delta_t;
                        dX_d  = (X_d  - [simul_X_d_x(n-1); simul_X_d_y(n-1)])./delta_t;
                        ddX_d = (dX_d  - [simul_dX_d_x(n-1); simul_dX_d_y(n-1)])./delta_t;
                    end
                end
                % Get Dynamics
                    J   = get_Jacobian(q(1),q(2));
                    dJ  = get_dJacobian(q(1),q(2));     % (J-pre_J)/delta_t
                    X   = get_Kinematics(q(1),q(2));
                    dX  = J*dq;                         % (X-pre_X)/delta_t
                    D   = [m2*L2^2 + 2*m2*L1*L2*cos(q(2)) + L1^2*(m1+m2), m2*L2^2+m2*L1*L2*cos(q(2));
                           m2*L2^2 + m2*L1*L2*cos(q(2)),                  m2*L2^2];
                    H   = [-m2*L1*L2*dq(2)^2*sin(q(2)) - 2*m2*L1*L2*dq(1)*dq(2)*sin(q(2));
                           m2*L1*L2*dq(1)^2*sin(q(2))];
                    C   = [m2*L2*g*cos(q(1)+q(2)) + L1*g*cos(q(1))*(m1+m2);
                           m2*g*L2*cos(q(1)+q(2))];
                % Controller Dynamics
                    error       = X_d - X;
                    errorSum    = errorSum + error*delta_t;
                    u           = ddX_d + Kv(i)*(dX_d - dX) + Kp(i)*(X_d - X) + Ki(i)*errorSum;
                    
                    ddq_ref  = inv(J)*(u - dJ*dq);
                    
                    tq_ctrl  = D*ddq_ref + H + C*0.8;  % 0.8 : 중력 보상 오차
                % Robot Model
                    % Inverse Dynamics
                        tq      = tq_ctrl;
                        tq1     = tq(1);    tq2 = tq(2);
                        [t, y]  = ode45('two_link' , [0 delta_t], [q(1); dq(1); q(2); dq(2)]);
                        index   = length(y);
                        q       = [y(index,1); y(index,3)];
                        dq      = [y(index,2); y(index,4)];
                        
                % Save Data
                    simul_time(n)   = time;       % [sec]
                    simul_q1(n)    = q(1);      % [rad]
                    simul_q2(n)    = q(2);      % [rad]
                    simul_dq1(n)   = dq(1);     % [rad/s]
                    simul_dq2(n)   = dq(2);     % [rad/s]
                    simul_X_x(n)   = X(1);      % [m]
                    simul_X_y(n)   = X(2);      % [m]
                    simul_dX_x(n)  = dX(1);     % [m/s]
                    simul_dX_y(n)  = dX(2);     % [m/s]   
                    simul_X_d_x(n) = X_d(1);    % [m]
                    simul_X_d_y(n) = X_d(2);    % [m]
                    simul_dX_d_x(n)= dX_d(1);   % [m/s]
                    simul_dX_d_y(n)= dX_d(2);   % [m/s]  
                    n             = n + 1;
            end
    end
%
%% Simulation Result Graph
    if(flag_Draw == 1)
        font_size_label         = 20;
        font_size_title         = 25;
        linewidth_current       = 3;
        linewidth_target        = 5;
        
        if(flag_Draw_Robot == 1)
            % Draw Robot
                x1  = L1*cos(init_q1);          % [m]   : Joint 1 X-axis Position
                y1  = L1*sin(init_q1);          % [m]   : Joint 1 Y-axis Position
                x2  = L2*cos(init_q1+init_q2);  % [m]   : Joint 2 X-axis Position
                y2  = L2*sin(init_q1+init_q2);  % [m]   : Joint 2 Y-axis Position

                
                FG1 = figure('Position', [200 300 500 450], 'Color', [1 1 1]);
                    AX = axes('parent', FG1); hold on
                    
                    Px1 = [0    x1];    Py1 = [0   y1];
                    Px2 = [x1   x1+x2]; Py2 = [y1  y1+y2];

                    p1 = plot(Px1, Py1, '-ob', 'LineWidth', linewidth_current);
                    p2 = plot(Px2, Py2, '-or', 'LineWidth', linewidth_current);

                    axis([-1.0 1.0 -1.6 0.4]);
                    grid on;
                xlabel('X-axis (m)', 'fontsize', font_size_label)
                ylabel('y-axis (m)', 'fontsize', font_size_label) 
                title('2-DOF Robot', 'fontsize', font_size_title)
                
                n =1;
                for (time = start_t: delta_t: finish_t)
                    q1  = simul_q1(n);  q2 = simul_q2(n);
                    
                    x1  = L1*cos(q1);         
                    y1  = L1*sin(q1);          
                	x2  = L2*cos(q1+q2);  
                    y2  = L2*sin(q1+q2); 
                    
                    Px1 = [0    x1];    Py1 = [0   y1];
                    Px2 = [x1   x1+x2]; Py2 = [y1  y1+y2];
                    set(p1,'XData',Px1,'YData',Py1)
                    set(p2,'XData',Px2,'YData',Py2)
                    drawnow
                    n = n+1;
                    
                    cmd = sprintf("Time : %2.2f",time);
                    clc
                    disp(cmd)
                end
        end
        if(flag_Draw_Graph == 1)
            % Draw Angle
                FG2 = figure('Position', [800 700 600 300], 'Color', [1 1 1]);
                    plot(simul_time,simul_X_d_x, ':r', 'linewidth', linewidth_target); hold on;
                    plot(simul_time,simul_X_d_y, ':b', 'linewidth', linewidth_target); hold on;
                   
                    plot(simul_time,simul_X_x, 'r', 'linewidth', linewidth_current); hold on;
                    plot(simul_time,simul_X_y, 'b', 'linewidth', linewidth_current); hold on;
                    
                    axis([start_t finish_t -1.25 1]);
                    xticks([start_t:1:finish_t])
                    yticks([-1:0.25:1])
                    grid on
                    
                    legend({'tar_x','tar_y', 'cur_x', 'cur_y'}, 'location', 'best', 'orientation', 'horizontal', 'fontsize', 10);
                xlabel('time (s)', 'fontsize', font_size_label)
                ylabel('Position (m)', 'fontsize', font_size_label) 
                title('Cartesian Space PID CTM Controller', 'fontsize', font_size_label)
            % Draw Velocity
                  FG3 = figure('Position', [800 100 600 300], 'Color', [1 1 1]);
                    plot(simul_time,simul_dX_d_x, ':r', 'linewidth', linewidth_target); hold on;
                    plot(simul_time,simul_dX_d_y, ':b', 'linewidth', linewidth_target); hold on;
                    
                    plot(simul_time,simul_dX_x, 'r', 'linewidth', linewidth_current); hold on;
                    plot(simul_time,simul_dX_y, 'b', 'linewidth', linewidth_current); hold on;
                    
                    axis([start_t finish_t -1.25 1.25]);
                    xticks([start_t:1:finish_t])
                    yticks([-1.25:0.25:1.25])
                    grid on
                    
                    legend({'tar_{dx}','tar_{dy}', 'cur_{dx}', 'cur_{dy}'}, 'location', 'best', 'orientation', 'horizontal', 'fontsize', 10);
                xlabel('time (s)', 'fontsize', font_size_label)
                ylabel('Velocity (m/s)', 'fontsize', font_size_label) 
                title('Cartesian Space PID CTM Controller', 'fontsize', font_size_label)
        end
    end                                           