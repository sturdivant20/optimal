%% Optimal HW4 - Problem 1 | Daniel Sturdivant
clc; close all; clear;
fprintf("<strong>PROBLEM 1</strong>\n");
f = figure;
tbs = uitabgroup(Parent=f);
tab(1) = uitab(Parent=tbs, Title="Sys");
tab(2) = uitab(Parent=tbs, Title="EKF");
tab(3) = uitab(Parent=tbs, Title="UKF");

% parameters
J = 2.5;
m = 1.6;
l = 1;
g = 9.80;
b = 1.25;
T = 12;

% freq and time
f = 100;
t = 0:1/f:10;
len = length(t);


%% PART A
fprintf("<strong>(a)</strong>\n");

theta = zeros(1,len);
omega = zeros(1,len);
alpha = zeros(1,len);
y = zeros(1,len);
y(1) = deg2rad(randn);

for k = 2:len
    
    % inputs
    dt = t(k) - t(k-1);
    F = 5 + sqrt(2)*randn;

    % system
    alpha(k) = 1/J * (T + F*l*cos(theta(k-1)) - b*omega(k-1)^3 - m*g*l*sin(theta(k-1)));
    omega(k) = omega(k-1) + alpha(k-1)*dt;
    theta(k) = wrapToPi(theta(k-1) + omega(k)*dt + 0.5*alpha(k)*dt^2);

    % measurements
    y(k) = theta(k) + deg2rad(randn);

end

tl = tiledlayout(2,1,Parent=tab(1), TileSpacing="tight");
axes(Parent=tl);
hold on;
plot(t, y, '.', MarkerSize=10);
plot(t, theta, LineWidth=3);
legend("Measurement", "Theta");
xlabel("Time [s]");
ylabel("Angle [rad]")

nexttile;
hold on;
plot(t, omega);
legend("Omega");
xlabel("Time [s]");
ylabel("Anglular Velocity [rad/s]");
ylim([-1.5,2]);


%% PART B
fprintf("<strong>(b)</strong>\n");

% x = [theta; omega; F]
alpha_ekf = zeros(1,len);
x_ekf = zeros(3,len);
P = zeros(3,3,len);
Q = diag([1,1,1].*1e-6);
R = 1;

C = [1, 0, 0];
P(:,:,1) = eye(3);
x_ekf(3,1) = 5;
L = P(:,:,1)*C'*(C*P(:,:,1)*C' + R)^-1;
P(:,:,1) = (eye(3) - L*C) * P(:,:,1);
x_ekf(:,1) = x_ekf(:,1) + L*(y(1) - C*x_ekf(:,1));


for k = 2:len

    dt = t(k) - t(k-1);

    % predict
    A = [                                                              0,                     1,   0; ...
         -(l*x_ekf(3,k-1)*sin(x_ekf(1,k-1)) + m*g*l*cos(x_ekf(1,k-1)))/J, -3*b*x_ekf(2,k-1)^2/J, l*cos(x_ekf(1,k-1))/J; ...
                                                                       0,                     0,   0];
%     A = expm(A*dt);
    A = eye(3) + A*dt;

    alpha_ekf(k) = 1/J * (T + x_ekf(3,k-1)*l*cos(x_ekf(1,k-1)) - b*x_ekf(2,k-1)^3 - m*g*l*sin(x_ekf(1,k-1)));
    x_ekf(1,k) = wrapToPi(x_ekf(1,k-1) + x_ekf(2,k-1)*dt + 0.5*alpha_ekf(k)*dt^2);
    x_ekf(2,k) = x_ekf(2,k-1) + alpha_ekf(k)*dt;
    x_ekf(3,k) = x_ekf(3,k-1);
    P(:,:,k) = A*P(:,:,k-1)*A' + Q;

    % correct
    C = [1, 0, 0];
    L = P(:,:,k)*C'*(C*P(:,:,k)*C' + R)^-1;
    P(:,:,k) = (eye(3) - L*C) * P(:,:,k);
    x_ekf(:,k) = x_ekf(:,k) + L*(y(k) - C*x_ekf(:,k));

end

tl = tiledlayout(2,1,Parent=tab(2), TileSpacing="tight");
axes(Parent=tl);
hold on;
plot(t, y, '.', MarkerSize=10);
plot(t, theta, LineWidth=3);
plot(t, x_ekf(1,:), LineWidth=1.5);
legend("Measurement", "Theta", "EKF \theta");
xlabel("Time [s]");
ylabel("Angle [rad]")

nexttile;
hold on;
plot(t, omega);
plot(t, x_ekf(2,:));
legend("Omega", "EKF \omega");
xlabel("Time [s]");
ylabel("Anglular Velocity [rad/s]");
ylim([-1.5,2]);


%% PART C
fprintf("<strong>(c)</strong>\n");

% x = [theta; omega; F]
n = 3;
alpha = 1e-3;
beta = 2;
kappa = 0;
lambda = alpha^2 * (n + kappa) - n;

alpha_ukf = zeros(2*n+1,len);
x_ukf = zeros(n,len);
y_ukf = zeros(1,len);
P_ukf = zeros(n,n,len);

Q = diag([1,1,1].*1e-6);
R = 1;
C = zeros(1,n);
C(1) = 1;

xHat = zeros(n, 2*n + 1);   % sigma point state mean
yHat = zeros(1, 2*n + 1);   % sigma point measurement mean
sig = zeros(n, 2*n + 1);    % sigma points
Wm = zeros(1, 2*n + 1);     % sigma point measurement weights
Wc = zeros(1, 2*n + 1);     % sigma point covariance weights

% first update
P_ukf(:,:,1) = eye(n)+Q;
x_ukf(3,1) = 5;
Pxy = zeros(n,1);
Pyy = 1+R;
L = Pxy * Pyy^-1;
P_ukf(:,:,1) = P_ukf(:,:,1) - L*Pyy*L';
x_ukf(:,1) = x_ukf(:,1) + L*(y(1) - y_ukf(1));

for k = 2:len

    % k is index of system
    % i is index of sigma point 
    
    dt = t(k) - t(k-1);

    % square root of covariance
    P_sqrt = chol((n + lambda) * P_ukf(:,:,k-1));
    Pxy = zeros(n,1);
    Pyy = 0;

    % first sigma point
    sig(:,1) = x_ukf(:,k-1);
    Wm(1) = lambda / (n + lambda);
    Wc(1) = Wm(1) + (1 - alpha^2 + beta);

    alpha_ukf(1,k) = 1/J * (T + sig(3,1)*l*cos(sig(1,1)) - b*sig(2,1)^3 - m*g*l*sin(sig(1,1)));
    xHat(1,1) = wrapToPi(sig(1,1) + sig(2,1)*dt + 0.5*alpha_ukf(1,k)*dt^2);
    xHat(2,1) = sig(2,1) + alpha_ukf(1,k) * dt;
    xHat(3,1) = sig(3,1);

    yHat(1) = C * xHat(:,1);
    x_ukf(:,k) = x_ukf(:,k) + Wm(1) * xHat(:,1);
    y_ukf(k) = y_ukf(k) + Wm(1) * yHat(1);

    for i = 1:2*n
        % sigma points
        if i < 4
            sig(:,i+1) = sig(:,1) + P_sqrt(:,i);
        else
            sig(:,i+1) = sig(:,1) - P_sqrt(:,i-n);
        end
        Wm(i+1) = 1 / (2*(n+lambda));
        Wc(i+1) = Wm(i+1);

        % propogated sigma point mean
        alpha_ukf(i+1,k) = 1/J * (T + sig(3,i+1)*l*cos(sig(1,i+1)) - b*sig(2,i+1)^3 - m*g*l*sin(sig(1,i+1)));
        xHat(1,i+1) = wrapToPi(sig(1,i+1) + sig(2,i+1)*dt + 0.5*alpha_ukf(i+1,k)*dt^2);
        xHat(2,i+1) = sig(2,i+1) + alpha_ukf(i+1,k) * dt;
        xHat(3,i+1) = sig(3,i+1);
        x_ukf(:,k) = x_ukf(:,k) + Wm(i+1) * xHat(:,i+1);

        % propagated sigma point measurement mean
        yHat(i+1) = C * xHat(:,i+1);
        y_ukf(k) = y_ukf(k) + Wm(i+1) * yHat(i+1);
    end

    % propogated system covariance
    for i = 1:(2*n + 1)
        P_ukf(:,:,k) = P_ukf(:,:,k) + Wc(i) * ( (xHat(:,i)-x_ukf(:,k)) * (xHat(:,i)-x_ukf(:,k))' );
        Pyy = Pyy + Wc(i) * ( (yHat(i)-y_ukf(k)) * (yHat(i)-y_ukf(k))' );
        Pxy = Pxy + Wc(i) * ( (xHat(:,i)-x_ukf(:,k)) * (yHat(i)-y_ukf(k))' );
    end
    Pyy = Pyy + R;
    P_ukf(:,:,k) = P_ukf(:,:,k) + Q;

    % correction
    L = Pxy * Pyy^-1;
    P_ukf(:,:,k) = P_ukf(:,:,k) - L*Pyy*L';
    x_ukf(:,k) = x_ukf(:,k) + L*(y(k) - y_ukf(k));

end

tl = tiledlayout(2,1,Parent=tab(3), TileSpacing="tight");
axes(Parent=tl);
hold on;
plot(t, y, '.', MarkerSize=10);
plot(t, theta, LineWidth=3);
plot(t, x_ukf(1,:), LineWidth=1.5);
legend("Measurement", "Theta", "UKF \theta");
xlabel("Time [s]");
ylabel("Angle [rad]")

nexttile;
hold on;
plot(t, omega);
plot(t, x_ukf(2,:));
legend("Omega", "UKF \omega");
xlabel("Time [s]");
ylabel("Anglular Velocity [rad/s]")
ylim([-1.5,2]);



