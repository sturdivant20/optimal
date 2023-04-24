%% Optimal HW4 - Problem 1 | Daniel Sturdivant
clc; close all; clear;
fprintf("<strong>PROBLEM 1</strong>\n");

%%

J = 2.5;
m = 1.6;
l = 1;
g = 9.80;
b = 1.25;
T = 12;

t = linspace(0,30,3000);
theta = zeros(1,3000);
omega = zeros(1,3000);
alpha = zeros(1,3000);
y = zeros(1,3000);
y(1) = deg2rad(randn);

for i = 2:3000
    
    % inputs
    dt = t(i) - t(i-1);
    F = 5 + sqrt(2)*randn;

    % system
    alpha(i) = 1/J * (T + F*l*cos(theta(i-1)) - b*omega(i-1)^3 - m*g*l*sin(theta(i-1)));
    omega(i) = omega(i-1) + alpha(i-1)*dt;
    theta(i) = wrapToPi(theta(i-1) + omega(i)*dt + 0.5*alpha(i-1)*dt^2);

    % measurements
    y(i) = theta(i) + deg2rad(randn);

end


% x = [theta; omega; F, b_theta, b_omega]
x = zeros(5,3000);
P = zeros(5,5,3000);
Q = diag([1, 1, 1, 1, 1]);
R = 1;

C = [1, 0, 0, 0, 0];
P(:,:,1) = eye(5);
x(3,1) = 5;
L = P(:,:,1)*C'*(C*P(:,:,1)*C' + R)^-1;
P(:,:,1) = (eye(5) - L*C) * P(:,:,1);
x(:,1) = x(:,1) + L*(y(1) - C*x(:,1));


for i = 2:3000

    dt = t(i) - t(i-1);

    % predict
    A = [1,   dt,                 0,               1,  dt; ...
         0,    1, l/J*cos(x(1,i-1)), -b/J*x(1,i-1)^3, -m*g*l*sin(x(1,i-1)); ...
         0,    0,                 1,               0,   0; ...
         0,    0,                 0,               1,  dt; ...
         0,    0,                 0,               0,   1];
    B = [0;
         1/J;
         0;
         0;
         0];

    x(:,i) = A*x(:,i-1); % + B*T;
    P(:,:,i) = A*P(:,:,i-1)*A' + Q;

    % correct
    C = [1, 0, 0, 0, 0];
    L = P(:,:,i)*C'*(C*P(:,:,i)*C' + R)^-1;
    P(:,:,i) = (eye(5) - L*C) * P(:,:,i);
    x(:,i) = x(:,i) + L*(y(i) - C*x(:,i));

end

figure;
hold on;
% plot(t, y, '.', MarkerSize=10);
plot(t, theta, LineWidth=3);
plot(t, x(1,:), LineWidth=1.5);
% legend("Measurement", "Theta", "KF \theta");
legend("Theta", "KF \theta");
xlabel("Time [s]");
ylabel("Angle [rad]")

figure
hold on;
plot(t, omega);
plot(t, x(2,:));
legend("Omega", "KF \omega");
xlabel("Time [s]");
ylabel("Anglular Velocity [rad/s]")

