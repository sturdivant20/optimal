%% Optimal HW3 - Problem 4 | Daniel Sturdivant
clc; clear; close all;
fprintf("<strong>PROBLEM 4</strong>\n");

% f = figure('units','normalized','position',[0.1 0.1 0.8 0.8]);
% f = figure(Units='normalized', Position=[1.1, 0.5, 0.8, 0.4]);
f = figure(Units='normalized', Position=[3.0, 0.5, 1.2, 0.4]);
tbs = uitabgroup(Parent=f);
tab(1) = uitab(Title="A", Parent=tbs);
tab(2) = uitab(Title="B", Parent=tbs);
tab(3) = uitab(Title="D", Parent=tbs);

dt = 0.1;
M = 5*10 + 1; % 5 min
t = linspace(0,5,M);


%% PART A
fprintf("\n<strong>(a)</strong>\n");

% continuous state
%  -> x_dot = Ax + Bu
A = [-2.62, 12; ...
     -0.96, -2];
B = [14; ...
      1];
C = [1, 0];

% Noise (guessing process noise)
R = 0.1^2;
Q = 0.05;

x = zeros(2,M);
y = zeros(1,M);
for k = 1:M
    if k == 1
        x(:,k) = [0;0]; % u=1
        y(k) = R*randn;
    else
        x_dot = A * x(:,k-1) + B * 1;
        x(:,k) = x(:,k-1) + x_dot*dt;
        y(k) = C * x(:,k) + R*randn;
    end
end

% bryson's rule
S1 = [        -A, B*Q*B'; ...
      zeros(2,2),     A'];
S2 = eye(4) + S1*dt;
% S2 = expm(S1);
Ad = S2(3:4,3:4)';
Bd = B*dt;
% Qd = Ad * S2(1:2,3:4);
Qd = [2, -1; -1, 0.5];

% Discrete measurement
%  -> y(k) = C*x(k) + v
C = [1, 0];

% simulation
x_hat = zeros(2,M);
P = zeros(2,2,M);
L = zeros(2,M);

for k = 1:M
    if k == 1
        % initialize kalman filter
        x_hat(:,k) = [0;0];
        P(:,:,k) = eye(2)*100;

    else
        % time update (propagate)
        x_hat(:,k) = Ad * x_hat(:,k-1); % only measure yaw rate
        P(:,:,k) = Ad * P(:,:,k-1) * Ad' + Qd;

    end

    % measurement update (correct)
    L(:,k) = P(:,:,k) * C' * (C * P(:,:,k) * C' + R)^-1;
    P(:,:,k) = (eye(2) - L(:,k) * C) * P(:,:,k);
    x_hat(:,k) = x_hat(:,k) + L(:,k) * (y(k) - C * x_hat(:,k));

end

poles = eig(Ad - L(:,end)*C);
fprintf("SS Poles = [%.3f%+.3fi; %.3f%+.3fi] \n", real(poles(1)), imag(poles(1)), real(poles(2)), imag(poles(2)));

tl = tiledlayout(2,1, Parent=tab(1), TileSpacing="tight");

axes(Parent=tl);
hold on;
plot(t, y, 'g.', MarkerSize=10, LineWidth=8);
plot(t, x(1,:), 'b--', LineWidth=2);
plot(t, x_hat(1,:), 'r', LineWidth=1.5);
grid on;
title("\bf{A) Yaw Rate}");

lg = legend("Meaurement", "State", "KF Estimate", Orientation="horizontal");
lg.Layout.Tile = "north";

nexttile;
hold on;
plot(t, x(2,:), 'b--', LineWidth=2);
plot(t, x_hat(2,:), 'r', LineWidth=1.5);
grid on;
title("\bf{A) Side Slip Angle}");
xlabel("Time [s]");

set(findall(gcf,'-property','FontSize'),'FontSize',16)
set(findall(gcf,'-property','Interpreter'),'Interpreter','latex')


%% PART B
fprintf("\n<strong>(b)</strong>\n");

% part a with different dynamic model but the same kalman filter
A2 = [-2.42,  4; ...
      -0.99, -2];
B2 = [18; ...
       1];

% Noise (guessing process noise)
R = 0.1^2;

x = zeros(2,M);
y = zeros(1,M);
for k = 1:M
    if k == 1
        x(:,k) = [0;0]; % u=1
        y(k) = R*randn;
    else
        x_dot = A2 * x(:,k-1) + B2 * 1;
        x(:,k) = x(:,k-1) + x_dot*dt;
        y(k) = C * x(:,k) + R*randn;
    end
end

% brysons rule
% S1 = [       -A2, B2*Q*B2'; ...
%       zeros(2,2),      A2'];
% S2 = eye(4) + S1*dt;
% Ad2 = S2(3:4,3:4)';
% Bd2 = B2*dt;
% Qd2 = Ad2 * S1(1:2,3:4);
Qd2 = [3, -1.2; -1.2, 0.5];

% simulation
x_hat = zeros(2,M);
P = zeros(2,2,M);
L = zeros(2,M);

for k = 1:M
    if k == 1

        % initialize kalman filter
        x_hat(:,k) = [0;0];
        P(:,:,k) = eye(2);

    else

        % time update (propagate) -> same as PART A
        x_hat(:,k) = Ad * x_hat(:,k-1); % only measure yaw rate
        P(:,:,k) = Ad * P(:,:,k-1) * Ad' + Qd2;

    end

    % measurement update (correct)
    L(:,k) = P(:,:,k) * C' * (C * P(:,:,k) * C' + R)^-1;
    P(:,:,k) = (eye(2) - L(:,k) * C) * P(:,:,k);
    x_hat(:,k) = x_hat(:,k) + L(:,k) * (y(k) - C * x_hat(:,k));

end

poles = eig(Ad - L(:,end)*C);
fprintf("SS Poles = [%.3f%+.3fi; %.3f%+.3fi] \n", real(poles(1)), imag(poles(1)), real(poles(2)), imag(poles(2)));

tl = tiledlayout(2,1, Parent=tab(2), TileSpacing="tight");

axes(Parent=tl);
hold on;
plot(t, y, 'g.', MarkerSize=10, LineWidth=8);
plot(t, x(1,:), 'b--', LineWidth=2);
plot(t, x_hat(1,:), 'r', LineWidth=1.5);
grid on;
title("\bf{B) Yaw Rate}");

lg = legend("Meaurement", "State", "KF Estimate", Orientation="horizontal");
lg.Layout.Tile = "north";

nexttile;
hold on;
plot(t, x(2,:), 'b--', LineWidth=2);
plot(t, x_hat(2,:), 'r', LineWidth=1.5);
grid on;
title("\bf{B) Side Slip Angle}");
xlabel("Time [s]");

set(findall(gcf,'-property','FontSize'),'FontSize',16)
set(findall(gcf,'-property','Interpreter'),'Interpreter','latex')


%% PART C
fprintf("\n<strong>(c)</strong>\n");

% now we have a measurement of slip angle
C3 = eye(2);
R3 = diag([0.1, 0.5].^2);
Q = 0.5;

fprintf("R = [%.3f, %.3f; %.3f, %.3f] \n", R3);


%% PART D
fprintf("\n<strong>(d)</strong>\n");

% brysons rule
S1 = [        -A, B*Q*B'; ...
      zeros(2,2),     A'];
% S2 = eye(4) + S1*dt;
S2 = expm(S1*dt);
Ad3 = S2(3:4,3:4)';
Bd3 = B*dt;
% Qd3 = Ad3 * S1(1:2,3:4);
Qd3 = [0.2, -0.01; -0.01, 0.02];
% Qd3 = [3, -1.2; -1.2, 0.5];

% simulation
x = zeros(2,M);
y = zeros(2,M);
x_hat = zeros(2,M);
P = zeros(2,2,M);
L = zeros(2,2,M);

for k = 1:M
    if k == 1
        % initialize system dynamics
        x(:,k) = [0;0]; % u=1
        y(:,k) = R*randn(2,1);

        % initialize kalman filter
        x_hat(:,k) = [0;0];
        P(:,:,k) = eye(2);

    else
        % system dynamics
        x(:,k) = Ad3 * x(:,k-1) + Bd3 * 1; % u=1
        y(:,k) = C3 * x(:,k) + R3*randn(2,1);

        % time update (propagate)
        x_hat(:,k) = Ad3 * x_hat(:,k-1); % only measure yaw rate
        P(:,:,k) = Ad3 * P(:,:,k-1) * Ad3' + Qd3;

    end

    % measurement update (correct)
    L(:,:,k) = P(:,:,k) * C3' * (C3 * P(:,:,k) * C3' + R3)^-1;
    P(:,:,k) = (eye(2) - L(:,:,k) * C3) * P(:,:,k);
    x_hat(:,k) = x_hat(:,k) + L(:,:,k) * (y(:,k) - C3 * x_hat(:,k));

end

poles = eig(Ad - L(:,:,end)*C3);
fprintf("SS Poles = [%.3f%+.3fi; %.3f%+.3fi] \n", real(poles(1)), imag(poles(1)), real(poles(2)), imag(poles(2)));

tl = tiledlayout(2,1, Parent=tab(3), TileSpacing="tight");

axes(Parent=tl);
hold on;
plot(t, y(1,:), 'g.', MarkerSize=10, LineWidth=8);
plot(t, x(1,:), 'b--', LineWidth=2);
plot(t, x_hat(1,:), 'r', LineWidth=1.5);
grid on;
title("\bf{D) Yaw Rate}");

nexttile;
hold on;
plot(t, y(2,:), 'g.', MarkerSize=10, LineWidth=8);
plot(t, x(2,:), 'b--', LineWidth=2);
plot(t, x_hat(2,:), 'r', LineWidth=1.5);
grid on;
title("\bf{D) Side Slip Angle}");
xlabel("Time [s]");

lg = legend("Measurement", "State", "KF Estimate", Orientation="horizontal");
lg.Layout.Tile = "north";
set(findall(gcf,'-property','FontSize'),'FontSize',16)
set(findall(gcf,'-property','Interpreter'),'Interpreter','latex')


%%

% exportgraphics(tab(1), "./media/p4_a.png");
% exportgraphics(tab(2), "./media/p4_b.png");
% exportgraphics(tab(3), "./media/p4_d.png");

