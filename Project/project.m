clc; clear; close all;

load("sim.mat");
mode = 'view';
% mode = 'mc';
saveFig = false;

n = 4;                      % number of states
m = size(y,1);              % number of measurements
s = size(x_sensor,1);       % number of sensors
L = size(t,2);              % number of time points
N = 1000;                   % number of particles
K = 1000;                   % number of monte carlo runs

R = diag([3,3,3,3]);    % measurement noise [deg]
RR = 3;
Q = diag([3, 3, 3, 3].*110);   % process noise [m, m, m/s, m/s]

% RUN SIMULATION
if strcmp(mode, 'view')
    % view plot gif of one run

    % INITIALIZE
    x = zeros(n,L);
    P = zeros(n,n,L);
    
    % FIRST UPDATE
    y_(:,1) = y(:,1);
    x(:,1) = [pos(:,1); vel(:,1)];
    p = diag([10,10,1,1])*randn(4,N) + x(:,1);
    P(:,:,1) = eye(4).*10;
    w = ones(1,N) ./ N;

    f = figure(Units='normalized', Position=[3.0, 0.4, 0.8, 0.4]);
    for i = 1:L
    
        if i > 1
        dt = t(i) - t(i-1);
            y_(:,i) = y(:,i) + 3*randn(m,1);
%             [x(1:2,i), R_ls] = ls(x_sensor, y_(:,i), x(1:2,i-1), RR);
%             [p, x(:,i), w] = pf(p, w, x_sensor, y_(:,i), Q, RR, dt, i);
%             [x(:,i), P(:,:,i)] = ekfL(x(:,i-1), P(:,:,i-1), x(1:2,i), Q, R(1:2,1:2), dt, i);
%             [x(:,i), P(:,:,i)] = ekfT(x(:,i-1), P(:,:,i-1), y_(:,i), x_sensor, Q, R, dt, i);
%             [x(:,i), P(:,:,i)] = ukf(x(:,i-1), P(:,:,i-1), y_(:,i), x_sensor, Q./1e1, R.*1, dt, i);
            [p, x(:,i), P(:,:,i)] = enkf(p, x(:,i-1), P(:,:,i-1), y_(:,i), x_sensor, Q./2e3, R.*1, dt, i);
        end
        
        % plot
        plot(x_sensor(:,1), x_sensor(:,2), 'bo');
        hold on;
%         plot(p(1,:), p(2,:), '.');
        plot(pos(1,1:i), pos(2,1:i), 'g*');
        plot(x(1,1:i), x(2,1:i), 'rx');
        for j = 1:length(y(:,i))
            plot([x_sensor(j,1), x_sensor(j,1) - 1000*cosd(y(j,i))], ...
                 [x_sensor(j,2), x_sensor(j,2) - 1000*sind(y(j,i))], 'k');
        end
        hold off;
        xlim([-120, 120]);
        ylim([-120, 120]);
%         legend('Sensor', 'Particles', 'Emitter', 'Estimate', Location='northoutside', Orientation='horizontal');
        legend('Sensor', 'Emitter', 'Estimate', Location='northoutside', Orientation='horizontal');
        title("LS");
        title("EKF");
        title("UKF");
        title("EnKF");
        grid on;
        pause(0.1);
    
    end
%     xlim([-120, 120]);
%     ylim([-120, 120]);

elseif strcmp(mode, 'mc')
    % run monte carlo

    % INITIALIZE
    x_ls = zeros(n/2,L, K);
    x_ekf = zeros(n,L, K);
    x_ukf = zeros(n,L, K);
    x_enkf = zeros(n,L, K);
    P_ls = zeros(n/2,n/2,L, K);
    P_ekf = zeros(n,n,L, K);
    P_ukf = zeros(n,n,L, K);
    P_enkf = zeros(n,n,L, K);
    diff_X_ls = zeros(n/2,L,K);
    diff_X_ekf = zeros(n,L,K);
    diff_X_ukf = zeros(n,L,K);
    diff_X_enkf = zeros(n,L,K);
    theor_x_std_ls = zeros(L,K);
    theor_x_std_ekf = zeros(L,K);
    theor_x_std_ukf = zeros(L,K);
    theor_x_std_enkf = zeros(L,K);
    theor_y_std_ls = zeros(L,K);
    theor_y_std_ekf = zeros(L,K);
    theor_y_std_ukf = zeros(L,K);
    theor_y_std_enkf = zeros(L,K);

    for k = 1:K

        % FIRST UPDATE
        y_(:,1) = y(:,1);
        x_ls(:,1,k) = [pos(:,1)];
        x_ekf(:,1,k) = [pos(:,1); vel(:,1)];
        x_ukf(:,1,k) = [pos(:,1); vel(:,1)];
        x_enkf(:,1,k) = [pos(:,1); vel(:,1)];
        p = diag([10,10,1,1])*randn(4,N) + x_enkf(:,1,k);
        P_ls(:,:,1,k) = eye(2).*10;
        P_ekf(:,:,1,k) = eye(4).*10;
        P_ukf(:,:,1,k) = eye(4).*10;
        P_enkf(:,:,1,k) = eye(4).*10;
        w = ones(1,N) ./ N;

        for i = 1:L

            if i > 1
                dt = t(i) - t(i-1);
                y_(:,i) = y(:,i) + 3*randn(m,1);
                [x_ls(:,i,k), P_ls(:,:,i,k)] = ls(x_sensor, y_(:,i), x_ls(:,i-1,k), RR);
                [x_ekf(:,i,k), P_ekf(:,:,i,k)] = ekfT(x_ekf(:,i-1,k), P_ekf(:,:,i-1,k), y_(:,i), x_sensor, Q, R, dt, i);
                [x_ukf(:,i,k), P_ukf(:,:,i,k)] = ukf(x_ukf(:,i-1,k), P_ukf(:,:,i-1,k), y_(:,i), x_sensor, Q./1e1, R./1, dt, i);
                [p, x_enkf(:,i,k), P_enkf(:,:,i,k)] = enkf(p, x_enkf(:,i-1,k), P_enkf(:,:,i-1,k), y_(:,i), x_sensor, Q./2e3, R.*1, dt, i);
            end
            
        end % end run

        diff_X_ls(:,:,k) = pos - x_ls(:,:,k);
        diff_X_ekf(:,:,k) = [pos; vel] - x_ekf(:,:,k);
        diff_X_ukf(:,:,k) = [pos; vel] - x_ekf(:,:,k);
        diff_X_enkf(:,:,k) = [pos; vel] - x_ekf(:,:,k);
        theor_x_std_ls(:,k) = sqrt(P_ls(1,1,:,k));
        theor_x_std_ekf(:,k) = sqrt(P_ekf(1,1,:,k));
        theor_x_std_ukf(:,k) = sqrt(P_ukf(1,1,:,k));
        theor_x_std_enkf(:,k) = sqrt(P_enkf(1,1,:,k));
        theor_y_std_ls(:,k) = sqrt(P_ls(2,2,:,k));
        theor_y_std_ekf(:,k) = sqrt(P_ekf(2,2,:,k));
        theor_y_std_ukf(:,k) = sqrt(P_ukf(2,2,:,k));
        theor_y_std_enkf(:,k) = sqrt(P_enkf(2,2,:,k));

    end % end monte carlo

    exp_x_std_ls = std(diff_X_ls(1,:,:), [], 3);
    exp_x_std_ekf = std(diff_X_ekf(1,:,:), [], 3);
    exp_x_std_ukf = std(diff_X_ukf(1,:,:), [], 3);
    exp_x_std_enkf = std(diff_X_enkf(1,:,:), [], 3);
    exp_y_std_ls = std(diff_X_ls(2,:,:), [], 3);
    exp_y_std_ekf = std(diff_X_ekf(2,:,:), [], 3);
    exp_y_std_ukf = std(diff_X_ukf(2,:,:), [], 3);
    exp_y_std_enkf = std(diff_X_enkf(2,:,:), [], 3);
    exp_x_mean_ls = mean(diff_X_ls(1,:,:), 3);
    exp_x_mean_ekf = mean(diff_X_ekf(1,:,:), 3);
    exp_x_mean_ukf = mean(diff_X_ukf(1,:,:), 3);
    exp_x_mean_enkf = mean(diff_X_enkf(1,:,:), 3);
    exp_y_mean_ls = mean(diff_X_ls(2,:,:), 3);
    exp_y_mean_ekf = mean(diff_X_ekf(2,:,:), 3);
    exp_y_mean_ukf = mean(diff_X_ukf(2,:,:), 3);
    exp_y_mean_enkf = mean(diff_X_enkf(2,:,:), 3);
    theor_x_std_ls = mean(theor_x_std_ls,2)';
    theor_x_std_ekf = mean(theor_x_std_ekf,2)';
    theor_x_std_ukf = mean(theor_x_std_ukf,2)';
    theor_x_std_enkf = mean(theor_x_std_enkf,2)';
    theor_y_std_ls = mean(theor_y_std_ls,2)';
    theor_y_std_ekf = mean(theor_y_std_ekf,2)';
    theor_y_std_ukf = mean(theor_y_std_ukf,2)';
    theor_y_std_enkf = mean(theor_y_std_enkf,2)';
    theor_mean = zeros(1,L);

    f = figure(Units='normalized', Position=[3.0, 0.4, 1.2, 0.4]);
    tbs = uitabgroup(Parent=f);
    tab(1) = uitab(Parent=tbs, Title='LS Error');
    tab(2) = uitab(Parent=tbs, Title='EKF Error');
    tab(3) = uitab(Parent=tbs, Title='UKF Error');
    tab(4) = uitab(Parent=tbs, Title='EnKF Error');

    % ls plot
    tl = tiledlayout(2,1, Parent=tab(1), TileSpacing="tight");
    axes(Parent=tl);
    hold on;
    plot(t, theor_mean, 'k--', LineWidth=2);
    plot(t, exp_x_mean_ls, 'g--', LineWidth=2);
    plot(t, [theor_x_std_ls; -theor_x_std_ls], 'r', LineWidth=2);
    plot(t, [3.*exp_x_std_ls; -3.*exp_x_std_ls], 'b', LineWidth=2);
    grid on;
    title("\textbf{X Position}");
    ylabel("Error [m]");
    legend("$\mu_{theoretical}$", "$\mu_{emperical}$", "$\sigma_{theoretical}$", "", "$\sigma_{emperical}$", "", Location="northoutside", Orientation="horizontal");
    nexttile;
    hold on;
    plot(t, theor_mean, 'k--', LineWidth=2);
    plot(t, exp_y_mean_ls, 'g--', LineWidth=2);
    plot(t, [theor_y_std_ls; -theor_y_std_ls], 'r', LineWidth=2);
    plot(t, [3.*exp_y_std_ls; -3.*exp_y_std_ls], 'b', LineWidth=2);
    grid on;
    title("\textbf{Y Position}");
    ylabel("Error [m]");
    xlabel("Time [s]");
    set(findall(gcf,'-property','FontSize'),'FontSize',16)
    set(findall(gcf,'-property','Interpreter'),'Interpreter','latex')

    % ekf plot
    tl = tiledlayout(2,1, Parent=tab(2), TileSpacing="tight");
    axes(Parent=tl);
    hold on;
    plot(t, theor_mean, 'k--', LineWidth=2);
    plot(t, exp_x_mean_ekf, 'g--', LineWidth=2);
    plot(t, [theor_x_std_ekf; -theor_x_std_ekf], 'r', LineWidth=2);
    plot(t, [3.*exp_x_std_ekf; -3.*exp_x_std_ekf], 'b', LineWidth=2);
    grid on;
    title("\textbf{X Position}");
    ylabel("Error [m]");
    legend("$\mu_{theoretical}$", "$\mu_{emperical}$", "$\sigma_{theoretical}$", "", "$\sigma_{emperical}$", "", Location="northoutside", Orientation="horizontal");
    nexttile;
    hold on;
    plot(t, theor_mean, 'k--', LineWidth=2);
    plot(t, exp_y_mean_ekf, 'g--', LineWidth=2);
    plot(t, [theor_y_std_ekf; -theor_y_std_ekf], 'r', LineWidth=2);
    plot(t, [3.*exp_y_std_ekf; -3.*exp_y_std_ekf], 'b', LineWidth=2);
    grid on;
    title("\textbf{Y Position}");
    ylabel("Error [m]");
    xlabel("Time [s]");
    set(findall(gcf,'-property','FontSize'),'FontSize',16)
    set(findall(gcf,'-property','Interpreter'),'Interpreter','latex')

    % ukf plot
    tl = tiledlayout(2,1, Parent=tab(3), TileSpacing="tight");
    axes(Parent=tl);
    hold on;
    plot(t, theor_mean, 'k--', LineWidth=2);
    plot(t, exp_x_mean_ukf, 'g--', LineWidth=2);
    plot(t, [theor_x_std_ukf; -theor_x_std_ukf], 'r', LineWidth=2);
    plot(t, [3.*exp_x_std_ukf; -3.*exp_x_std_ukf], 'b', LineWidth=2);
    grid on;
    title("\textbf{X Position}");
    ylabel("Error [m]");
    legend("$\mu_{theoretical}$", "$\mu_{emperical}$", "$\sigma_{theoretical}$", "", "$\sigma_{emperical}$", "", Location="northoutside", Orientation="horizontal");
    nexttile;
    hold on;
    plot(t, theor_mean, 'k--', LineWidth=2);
    plot(t, exp_y_mean_ukf, 'g--', LineWidth=2);
    plot(t, [theor_y_std_ukf; -theor_y_std_ukf], 'r', LineWidth=2);
    plot(t, [3.*exp_y_std_ukf; -3.*exp_y_std_ukf], 'b', LineWidth=2);
    grid on;
    title("\textbf{Y Position}");
    ylabel("Error [m]");
    xlabel("Time [s]");
    set(findall(gcf,'-property','FontSize'),'FontSize',16)
    set(findall(gcf,'-property','Interpreter'),'Interpreter','latex')

    % enkf plot
    tl = tiledlayout(2,1, Parent=tab(4), TileSpacing="tight");
    axes(Parent=tl);
    hold on;
    plot(t, theor_mean, 'k--', LineWidth=2);
    plot(t, exp_x_mean_enkf, 'g--', LineWidth=2);
    plot(t, [theor_x_std_enkf; -theor_x_std_enkf], 'r', LineWidth=2);
    plot(t, [3.*exp_x_std_enkf; -3.*exp_x_std_enkf], 'b', LineWidth=2);
    grid on;
    title("\textbf{X Position}");
    ylabel("Error [m]");
    legend("$\mu_{theoretical}$", "$\mu_{emperical}$", "$\sigma_{theoretical}$", "", "$\sigma_{emperical}$", "", Location="northoutside", Orientation="horizontal");
    nexttile;
    hold on;
    plot(t, theor_mean, 'k--', LineWidth=2);
    plot(t, exp_y_mean_enkf, 'g--', LineWidth=2);
    plot(t, [theor_y_std_enkf; -theor_y_std_enkf], 'r', LineWidth=2);
    plot(t, [3.*exp_y_std_enkf; -3.*exp_y_std_enkf], 'b', LineWidth=2);
    grid on;
    title("\textbf{Y Position}");
    ylabel("Error [m]");
    xlabel("Time [s]");
    set(findall(gcf,'-property','FontSize'),'FontSize',16)
    set(findall(gcf,'-property','Interpreter'),'Interpreter','latex')

end


%% 
if saveFig
    exportgraphics(tab(1), "./media/ls_error.png")
    exportgraphics(tab(2), "./media/ekf_error.png")
    exportgraphics(tab(3), "./media/ukf_error.png")
    exportgraphics(tab(4), "./media/enkf_error.png")
end
