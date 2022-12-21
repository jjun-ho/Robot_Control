%% Termproject 03

clc,clear
close all
%% Set Simulation Parameters
    % Draw flag
        flag_Simul      = 1;

        flag_Draw       = 1;
        flag_Draw_Robot = 1;
        flag_Draw_Graph = 1;

    % Global Variable
        global Iz1 L1 g m1 r1 tq1

    % Simulation Parameters
        delta_t     = 0.005;        % [sec]     : Sampling Time - 200Hz
        start_t     = 0.000;        % [sec]     : Start Time
        finish_t    = 5.000;        % [sec]     : End Time 

        g           = 9.8148;       % [m/s^2]   : Gravitational Acceleration
    % Robot Parameters
        m1          = 0.2;                % [kg]      : Link Mass
        L1          = 0.5;                % [m]       : Link Length
        Iz1         = 0.05;         
        r1          = 0.1;      
        tq1        = 0.000;              % [Nm]      : Control Torque

        init_q1     = 0.00;               % [rad]     : Init Joint Angle
        init_dq1    = 0.00;            % [rad/s]   : Init Angular Velocity
        q1          = init_q1;         % [rad]     : Joint Angle
        dq1         = init_dq1;       % [rad/s]   : Angular Velocity
        final_q1    = 90*pi/180;      % [rad]     : Final Joint Angle
        final_dq1   = 30*pi/180;     % [rad/s]   : Final Angular Velocity

        init_q      = [  init_q1];
        init_dq     = [ init_dq1];
        q           = [       q1];
        dq          = [      dq1];
        final_q     = [ final_q1];
        final_dq    = [final_dq1];

    % Target Pos/Vel/Acc Parameters
        q1_d        = init_q1;       % [rad]     : Target Joint Angle
        dq1_d       = init_dq1;      % [rad/s]   : Target Angular Velocity
        ddq1_d      = 0;             % [rad/s^2] : Target Angular Acceleration

        q_d         = [q1_d];
        dq_d        = [dq1_d];
        ddq_d       = [ddq1_d];

    % Controller Gain
        Wn          = 20;                % [rad/s]   : natural frequency
        Kp          = [Wn^2];            % [Nm/rad]  : Propotional Gain
        zeta        = 1;                 % Damping ratio(감쇠비)
        Kv          = [2*zeta*Wn];       % [Nm*s/rad]: Derivative Gain
        Ki          = [400];
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
                for (i = 1:1)
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
                    G           = get_Gravity(q);
                % Controller
                    error       = q_d - q;
                    errorSum    = errorSum + error*delta_t;
                    u           = ddq_d + Kv(i)*(dq_d - dq) + Kp(i)*(q_d - q) + Ki(i)*errorSum;
                    I            = (m1*L1^2)/3;    % [kgm^2]   : Link Inertia
                    tq_ctrl    = I*u + G * 0.8;

                end
                

                % Robot Model  
                    % Inverse Dynamics
                        tq1      = tq_ctrl(1);
                        [t,y]   = ode45('one_link',[0 delta_t],[q(1); dq(1)]);
                        index   = length(y);
                        q(1)    = y(index,1);
                        dq(1)   = y(index,2);

                    % Save Data
                        simul_time(n)   = time;     % [sec]
                        simul_q(:,n)    = q;        % [rad]
                        simul_dq(:,n)   = dq;       % [rad/s]
                        simul_q_d(:,n)  = q_d;      % [rad]
                        simul_dq_d(:,n) = dq_d;     % [rad/s]
                        n               = n + 1;

                    % Process
                        cmd = sprintf("Loading...(%2.0f%%)", time/finish_t*100);
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
                init_x = L1*cos(init_q);
                init_y = L1*sin(init_q);
                
                FG1 = figure('Position',[200 300 500 450],'Color',[1 1 1]);
                    AX = axes('parent',FG1); hold on
                    
                    p = plot([0 0],[init_x init_y],'-ob','Linewidth',linewidth_current);
                    
                    axis([-1.5 1.5 -1.5 1.5]);
                    grid on
               xlabel('X-axis (m)', 'fontsize',font_size_label)
               ylabel('Y-axis (m)', 'fontsize',font_size_label)
               title('1-DOF Robot','fontsize',font_size_title)
               
               n = 1;
               for(time = start_t:delta_t:finish_t)
                   q = simul_q(n);
                   x = L1*cos(q); y = L1*sin(q);
                   Px = [0,x];  Py = [0,y];
                   set(p,'XData',Px,'YData',Py)
                   drawnow
                   n = n + 1;
                   
                   cmd = sprintf("Time : %2.2f",time);
                   clc
                   disp(cmd)
               end
        end
        if(flag_Draw_Graph == 1)
            % Draw Angle
                FG2 = figure('Position',[800 700 600 300],'Color',[1 1 1]);
                    plot(simul_time,simul_q_d*180/pi,':k','Linewidth',linewidth_target); hold on
                    plot(simul_time,simul_q*180/pi,' r','Linewidth',linewidth_current); hold on
                    
                    legend('Desired','Current')
                    axis([start_t finish_t 0 120]);
                    xticks([start_t:1:finish_t])
                    yticks([0:45:90])
                    grid on
                xlabel('time (s)', 'fontsize',font_size_label)
                ylabel('Angle (deg)', 'fontsize',font_size_label)
                title('Joint Space PID CTM Controller','fontsize',font_size_title)
                
            % Draw Angular Velocity
                FG3 = figure('Position',[800 100 600 300],'Color',[1 1 1]);
                    plot(simul_time,simul_dq_d*180/pi,':k','Linewidth',linewidth_target); hold on
                    plot(simul_time,simul_dq*180/pi,' r','Linewidth',linewidth_current); hold on
                    
                    legend('Desired','Current')
                    axis([start_t finish_t -90 90]);
                    xticks([start_t:1:finish_t])
                    yticks([-60 0 45/2 60])
                    grid on
                xlabel('time (s)', 'fontsize',font_size_label)
                ylabel('Anglular Velocity (deg/s)', 'fontsize',font_size_label)
                title('Joint Space PID CTM Controller','fontsize',font_size_title)
        end
    end