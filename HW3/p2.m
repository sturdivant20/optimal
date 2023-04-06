%% Optimal HW3 - Problem 2 | Daniel Sturdivant
clc; clear; close all;
fprintf("<strong>PROBLEM 2</strong>\n");

% f = figure('units','normalized','position',[0.1 0.1 0.8 0.8]);
% f = figure(Units='normalized', Position=[1.1, 0.5, 0.8, 0.4]);
f = figure(Units='normalized', Position=[3.0, 0.5, 1.2, 0.4]);
tbs = uitabgroup(Parent=f);
tab(1) = uitab(Title="A", Parent=tbs);
tab(2) = uitab(Title="B", Parent=tbs);
tab(3) = uitab(Title="C", Parent=tbs);

% part 2 data
data = importdata("data\data_hw3_2.txt");
t = data(:,1);
y = data(:,2);
M = length(y);

% x_dot = 0;
A = 0;
Ad = 1;

% y = x + v;
C = 1;

% v ~ N(0,1=R)
R = 1;


%% PART A
fprintf("\n<strong>(a)</strong>\n");

% kalman filter with 0 process noise (Qd = 0)
Q = 0;

x = zeros(1,M);
P = zeros(1,M);
L = zeros(1,M);
I = eye(1);

for k = 1:M

    if k == 1
        % initialize
        x(:,k) = 0;
        P(:,:,k) = 1;

    else
        % time update (propogate)
        x(k) = Ad * x(k-1);
        P(k) = Ad * P(k-1) * Ad' + Q;
    end

    % measurement update (correct)
    L(k) = P(k) * C' * (C * P(k) * C' + R)^-1;
    P(k) = (I - L(k) * C) * P(k);
    x(k) = x(k) + L(k) * (y(k) - C * x(k));

end

% Lss
Lss = L(:,end);
fprintf("Lss = %.4f \n", Lss);

tl = tiledlayout(2,1, Parent=tab(1), TileSpacing="tight");

axes(parent=tl);
hold on;
plot(t,y, 'b.', LineWidth=4);
plot(t,x, 'r', LineWidth=2.5);
title("\bf{A) Simulated Measurements and KF Response}")
legend(["Measurement (y)", "KF Estimate ($\hat{x}$)"])
grid on;

nexttile;
plot(t,L, LineWidth=3);
grid on;
title("\bf{A) Kalman Gain}")
xlabel("Time [s]")

set(findall(gcf,'-property','FontSize'),'FontSize',16)
set(findall(gcf,'-property','Interpreter'),'Interpreter','latex')


%% PART B
fprintf("\n<strong>(b)</strong>\n");

% tuning Q
qRange = logspace(-2, -4, 4);

x = zeros(4, M);
P = zeros(4, M);
L = zeros(4, M);
Pss = zeros(1,4);
Lss = zeros(1,4);
I = eye(1);

for i = 1:4
    Q = qRange(i);

    for k = 1:M
    
        if k == 1
            % initialize
            x(i,k) = 0;
            P(i,k) = 1;
    
        else
            % time update (propogate)
            x(i,k) = Ad * x(i,k-1);
            P(i,k) = Ad * P(i,k-1) * Ad' + Q;

        end
    
        % measurement update (correct)
        L(i,k) = P(i,k) * C' * (C * P(i,k) * C' + R)^-1;
        P(i,k) = (I - L(i,k) * C) * P(i,k);
        x(i,k) = x(i,k) + L(i,k) * (y(k) - C * x(i,k));
    
    end

    Pss(i) = ((Ad * P(i,k) * Ad' + Q)^-1 + C * R^-1 * C')^-1;
    Lss(i) = Pss(i) * C' * R^-1;
end

tl = tiledlayout(2,1, Parent=tab(2), TileSpacing="tight");

axes(Parent=tl);
hold on;
plot(t,y, 'b.', LineWidth=5);
plot(t,x, LineWidth=3);
title("\bf{B) Simulated Measurements and KF Response}")
legend(["Measurement (y)", "Q=0.01", "Q=0.001", "Q=0.001", "Q=0.0001"], Location="southoutside", Orientation="horizontal");
grid on;

nexttile;
hold on;
plot(0,0);
plot(t,L, LineWidth=3);
grid on;
title("\bf{B) Kalman Gain}")
xlabel("Time [s]")

set(findall(gcf,'-property','FontSize'),'FontSize',16)
set(findall(gcf,'-property','Interpreter'),'Interpreter','latex')


%% PART C
fprintf("\n<strong>(c)</strong>\n");

% first order low pass filter 
%  -> H(z) = sqrt(Qd) / (z - (1-sqrt(Qd))
Q = 1e-2;
yf = filter(sqrt(Q), [1, -1+sqrt(Q)], y, y(1));

ax = axes(Parent=tab(3));
hold on;
plot(t,y, 'b.', LineWidth=4)
plot(t,x(1,:), 'r', LineWidth=4);
plot(t,yf, 'g--', LineWidth=3);
title("\bf{C) KF vs. LPF Response}")
legend(["Filterd Measurements", "Kalman Filter", "Low Pass Filter"]);
grid on;
set(findall(gcf,'-property','FontSize'),'FontSize',16)
set(findall(gcf,'-property','Interpreter'),'Interpreter','latex')


%% PART D
fprintf("\n<strong>(d)</strong>\n");

% why be same?? LPF = KF?


%% 

% exportgraphics(tab(1), "./media/p2_a.png");
% exportgraphics(tab(2), "./media/p2_b.png");
% exportgraphics(tab(3), "./media/p2_c.png");

