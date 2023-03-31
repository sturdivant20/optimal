%% Optimal HW3 - Problem 3 | Danel Sturdivant
clc; clear; close all;

% f = figure('units','normalized','position',[0.1 0.1 0.8 0.8]);
% f = figure(Units='normalized', Position=[1.1, 0.5, 0.8, 0.4]);
f = figure(Units='normalized', Position=[3.0, 0.5, 1.2, 0.4]);
tbs = uitabgroup(Parent=f);
tab(1) = uitab(Title="A", Parent=tbs);
tab(2) = uitab(Title="A2", Parent=tbs);
tab(3) = uitab(Title="B", Parent=tbs);
tab(4) = uitab(Title="B2", Parent=tbs);
tab(5) = uitab(Title="C", Parent=tbs);
tab(6) = uitab(Title="C2", Parent=tbs);

% part 2 data
data = importdata("data\data_hw3_3.txt");
t = data.data(:,1);
east = data.data(:,2);
north = data.data(:,3);
psi = wrapTo2Pi(data.data(:,4));
gyro = data.data(:,5);
radar = data.data(:,6);
N = length(t);

clearvars data;

%% PART A
fprintf("\n<strong>(a)</strong>\n");

% kalman filter ?
I = eye(5);
x = zeros(5,N);
P = zeros(5,5,N);
L = zeros(5,5,N);
u = zeros(2,N);

R = eye(5);
Q = 5e-1 .* diag([1,5,10,0.01,0.01]);

for k = 1:N
    % control
    u(:,k) = [gyro(k); ...
          radar(k)];

    y = [east(k); ...
         north(k); ...
         psi(k); ...
         0; ...
         0];

    if k == 1
        % initialize
        x(:,1) = y;
        P(:,:,k) = 100.*eye(5);

    else
        % time update (propagation)
        dt = t(k) - t(k-1);

        A = [1, 0, 0,  0 , -dt * sin(psi(k)); ...
             0, 1, 0,  0 , -dt * cos(psi(k)); ...
             0, 0, 1, -dt,  0; ...
             0, 0, 0,  1 ,  0; ...
             0, 0, 0,  0 ,  1];

        B = [ 0, dt * sin(psi(k)); ...
              0, dt * cos(psi(k)); ...
             dt, 0; ...
              0, 0; ...
              0, 0];

        x(:,k) = A * x(:,k-1) + B * u(:,k);
        P(:,:,k) = A * P(:,:,k-1) * A' + Q;

    end

    % measurement update (correct)
    C = eye(5);

    L(:,:,k) = P(:,:,k) * C' * (C * P(:,:,k) * C' + R)^-1;
    P(:,:,k) = (I - L(:,:,k) * C) * P(:,:,k);
    x(:,k) = x(:,k) + L(:,:,k) * (y - C * x(:,k));

end

ax = axes(Parent=tab(1));
hold on;
plot(east, north, 'b.', LineWidth=8, MarkerSize=10);
plot(x(1,:), x(2,:), 'r', LineWidth=2);
grid on;
legend("Measurement", "KF Estimate");
title("\bf{A) Trajectory Estimation}");
set(findall(gcf,'-property','FontSize'),'FontSize',16)
set(findall(gcf,'-property','Interpreter'),'Interpreter','latex')

tl = tiledlayout(3,2, Parent=tab(2), TileSpacing="tight");

axes(Parent=tl);
hold on;
plot(t, x(1,:), 'r', LineWidth=2);
plot(t, east, 'b--', LineWidth=1.5);
grid on;
ylabel("Position [m]");
title("\bf{East Position}");

nexttile;
hold on;
plot(t, x(4,:), 'r', LineWidth=2);
grid on;
ylabel("Bias [rad/s]");
title("\bf{Gyro Bias}");

nexttile;
hold on;
plot(t, x(2,:), 'r', LineWidth=2);
plot(t, north, 'b--', LineWidth=1.5);
grid on;
ylabel("Position [m]");
title("\bf{North Position}");

nexttile;
hold on;
plot(t, x(5,:), 'r', LineWidth=2);
grid on;
ylabel("Bias [m/s]");
title("\bf{Radar Bias}");

nexttile;
hold on;
plot(t, x(3,:), 'r', LineWidth=2);
plot(t, psi, 'b--', LineWidth=1.5);
grid on;
ylabel("Angle [rad]");
title("\bf{Heading Angle}");

lg = legend("KF Estimate", " Measurement", Location="northoutside", Orientation="horizontal");
lg.Layout.Tile = 'north';
set(findall(gcf,'-property','FontSize'),'FontSize',16)
set(findall(gcf,'-property','Interpreter'),'Interpreter','latex')


%% PART B
fprintf("\n<strong>(b)</strong>\n");

% recursive least squares ?

x1 = zeros(7,N);
P1 = zeros(7,7,N);

R = eye(5);
Q = 5e-1 .* diag([1,5,10,0.01,0.01]);

for k = 1:N
    dt = 0.2;

    % measurement matrix
    H = [1, 0, 0,  0, dt*sin(psi(k)),   0, -dt*sin(psi(k)); ...
         0, 1, 0,  0, dt*cos(psi(k)),   0, -dt*cos(psi(k)); ...
         0, 0, 1, dt,              0, -dt,               0; ...
         0, 0, 0,  1,              0,   0,               0; ...
         0, 0, 0,  0,              1,   0,               0];

    % measurement
    y = [east(k); ...
         north(k); ...
         psi(k); ...
         gyro(k); ...
         radar(k)];

    if k == 1
        % initialization
        G = 100.*eye(7);
        x1(:,k) = [y; ...
                   0; ...
                   0];
        P1(:,:,k) = G;

    else
        % recusive propagation
        G = G + H'*R^-1*H;
        dx = G^-1 * H' * R^-1 * (y - H * x1(:,k-1));
        x1(:,k) = x1(:,k-1) + dx;
        P1(:,:,k) = G;
    end
end


ax = axes(Parent=tab(3));
hold on;
plot(east, north, 'b.', LineWidth=8, MarkerSize=10);
plot(x1(1,:), x1(2,:), 'r', LineWidth=2);
grid on;
legend("Measurement", "LS Estimate");
title("\bf{A) Trajectory Estimation}");
set(findall(gcf,'-property','FontSize'),'FontSize',16)
set(findall(gcf,'-property','Interpreter'),'Interpreter','latex')

tl = tiledlayout(3,2, Parent=tab(4), TileSpacing="tight");

axes(Parent=tl);
hold on;
plot(t, x1(1,:), 'r', LineWidth=2);
plot(t, east, 'b--', LineWidth=1.5);
grid on;
ylabel("Position [m]");
title("\bf{East Position}");

nexttile;
hold on;
plot(t, x1(6,:), 'r', LineWidth=2);
grid on;
ylabel("Bias [rad/s]");
title("\bf{Gyro Bias}");

nexttile;
hold on;
plot(t, x1(2,:), 'r', LineWidth=2);
plot(t, north, 'b--', LineWidth=1.5);
grid on;
ylabel("Position [m]");
title("\bf{North Position}");

nexttile;
hold on;
plot(t, x1(7,:), 'r', LineWidth=2);
grid on;
ylabel("Bias [m/s]");
title("\bf{Radar Bias}");

nexttile;
hold on;
plot(t, x1(3,:), 'r', LineWidth=2);
plot(t, psi, 'b--', LineWidth=1.5);
grid on;
ylabel("Angle [rad]");
title("\bf{Heading Angle}");

lg = legend("LS Estiamte", "Measurement", Location="northoutside", Orientation="horizontal");
lg.Layout.Tile = 'north';
set(findall(gcf,'-property','FontSize'),'FontSize',16)
set(findall(gcf,'-property','Interpreter'),'Interpreter','latex')


%% PART C
fprintf("\n<strong>(c)</strong>\n");

I = eye(5);
x2 = zeros(5,N);
P2 = zeros(5,5,N);
L2 = zeros(5,5,N);
u2 = zeros(2,N);

R = eye(5);
Q = 5e-1 .* diag([1,5,10,0.01,0.01]);

for k = 1:N

    % system dynamics
    dt = 0.2; 

    A = [1, 0, 0,  0 , -dt * sin(psi(k)); ...
         0, 1, 0,  0 , -dt * cos(psi(k)); ...
         0, 0, 1, -dt,  0; ...
         0, 0, 0,  1 ,  0; ...
         0, 0, 0,  0 ,  1];

    B = [ 0, dt * sin(psi(k)); ...
          0, dt * cos(psi(k)); ...
         dt, 0; ...
          0, 0; ...
          0, 0];

    C = eye(5);

    % control
    u2(:,k) = [gyro(k); ...
               radar(k)];

    if t(k) < t(end)-40
    
        % measurement
        y = [east(k); ...
             north(k); ...
             psi(k); ...
             0; ...
             0];

        if k == 1
            % initialize
            x2(:,1) = y;
            P2(:,:,k) = 100.*eye(5);
    
        else
            % time update (propagation)   
            x2(:,k) = A * x2(:,k-1) + B * u2(:,k);
            P2(:,:,k) = A * P2(:,:,k-1) * A' + Q;
    
        end
    
        % measurement update (correct)
        L2(:,:,k) = P2(:,:,k) * C' * (C * P2(:,:,k) * C' + R)^-1;
        P2(:,:,k) = (I - L2(:,:,k) * C) * P2(:,:,k);
        x2(:,k) = x2(:,k) + L2(:,:,k) * (y - C * x2(:,k));
        x2(3,k) = wrapTo2Pi(x2(3,k));

    else
        % raw integration
        x2(:,k) = A*x2(:,k-1) + B*u2(:,k);
        x2(3,k) = wrapTo2Pi(x2(3,k));

    end

end


ax = axes(Parent=tab(5));
hold on;
plot(east, north, 'b.', LineWidth=8, MarkerSize=10);
plot(x2(1,:), x2(2,:), 'r', LineWidth=2);
grid on;
legend("Measurement", "LS Estimate");
title("\bf{A) Trajectory Estimation}");
set(findall(gcf,'-property','FontSize'),'FontSize',16)
set(findall(gcf,'-property','Interpreter'),'Interpreter','latex')

tl = tiledlayout(3,2, Parent=tab(6), TileSpacing="tight");

axes(Parent=tl);
hold on;
plot(t, x2(1,:), 'r', LineWidth=2);
plot(t, east, 'b--', LineWidth=1.5);
grid on;
ylabel("Position [m]");
title("\bf{East Position}");
xlim([250,300]);

nexttile;
hold on;
plot(t, x2(4,:), 'r', LineWidth=2);
grid on;
ylabel("Bias [rad/s]");
title("\bf{Gyro Bias}");
xlim([250,300]);

nexttile;
hold on;
plot(t, x2(2,:), 'r', LineWidth=2);
plot(t, north, 'b--', LineWidth=1.5);
grid on;
ylabel("Position [m]");
title("\bf{North Position}");
xlim([250,300]);

nexttile;
hold on;
plot(t, x2(5,:), 'r', LineWidth=2);
grid on;
ylabel("Bias [m/s]");
title("\bf{Radar Bias}");
xlim([250,300]);

nexttile;
hold on;
plot(t, x2(3,:), 'r', LineWidth=2);
plot(t, psi, 'b--', LineWidth=1.5);
grid on;
ylabel("Angle [rad]");
title("\bf{Heading Angle}");
xlim([250,300]);

lg = legend("KF Estimate", "Measurement", Location="northoutside", Orientation="horizontal");
lg.Layout.Tile = 'north';
set(findall(gcf,'-property','FontSize'),'FontSize',16)
set(findall(gcf,'-property','Interpreter'),'Interpreter','latex')


%%

% % exportgraphics(tab(1), "./media/p3_a.png");
% exportgraphics(tab(2), "./media/p3_a2.png");
% % exportgraphics(tab(3), "./media/p3_b.png");
% exportgraphics(tab(4), "./media/p3_b2.png");
% % exportgraphics(tab(5), "./media/p3_c.png");
% exportgraphics(tab(6), "./media/p3_c2.png");

