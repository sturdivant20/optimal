%% Optimal HW3 - Problem 1 | Daniel Sturdivant
clc; clear; close all;
fprintf("<strong>PROBLEM 1</strong>\n");

% f = figure('units','normalized','position',[0.1 0.1 0.8 0.8]);
% f = figure(Units='normalized', Position=[1.1, 0.5, 0.8, 0.4]);
f = figure(Units='normalized', Position=[3.0, 0.5, 1.2, 0.4]);
tbs = uitabgroup(Parent=f);
tab(1) = uitab(Title="A", Parent=tbs);
tab(2) = uitab(Title="C", Parent=tbs);
tab(3) = uitab(Title="C2", Parent=tbs);
tab(4) = uitab(Title="E", Parent=tbs);

% second order system with sample rate of 10Hz
dt = 0.1;

A = [ 0.0,  1.0; ...
     -1.0, -1.4];
Bw = [0.0; ...
      1.0];
C = [1.0, 0.0];

% noise characteristics
Qc = 2^2;
Rc = 1^2;


%% PART A
fprintf("\n<strong>(a)</strong>\n");

% simulate system for 100s
% x_dot = A*x + Bw*w
%     y = C*x + v

M = 1001;
x = zeros(2,M);
y = zeros(1,M);
x_dot = zeros(2,M);

for k = 1:M
    w = randn * 2;
    v = randn * 1;

    if k == 1
        % initialize
        x(:,k) = Bw*w;
        y(k) = v;
    else
        % measurement
        y(k) = C*x(:,k-1) + v;
    
        % dynamics
        x_dot = A*x(:,k-1) + Bw*w;
    
        % propogate
        x(:,k) = x(:,k-1) + dt*x_dot;
    end
    
end

% plot
t = linspace(0,(M-1)/10,M);

ax = axes(Parent=tab(1));
hold on
grid on
plot(t, y, '.b', LineWidth=4)
plot(t, x(1,:), 'r', LineWidth=3)
legend("Measurement (y)", "State ($x_1$)")
title("\bf{A) Simulated System Position State and Measurements}")
xlabel("Time [s]")
ylabel("Position")
set(findall(gcf,'-property','FontSize'),'FontSize',16)
set(findall(gcf,'-property','Interpreter'),'Interpreter','latex')

% figure
% plot(t, xp(2,:))


%% PART B
fprintf("\n<strong>(b)</strong>\n");

% constant dt therefore all discrete matrices are constant
I = eye(2);
% Ad = I + A*dt;
% Qd = Bw * Qc * Bw' * dt;
% Rd = Rc / dt;
% 
% fprintf("Ad = [%.3f %.3f; %.3f %.3f] \n", Ad)
% fprintf("Qd = [%.3f %.3f; %.3f %.3f] \n", Qd)
% fprintf("Rd = %.3f \n", Rd)

% BRYSON'S TRICK
S1 = [-A         , Bw*Qc*Bw'; ... 
       zeros(2,2), A'       ];
S2 = eye(4) + S1*dt;
% S2 = expm(S1*dt);
Ad = S2(3:4,3:4)';
Qd = Ad * S2(1:2,3:4);
Qd(1:3) = 0;
Rd = 1;

fprintf("Ad = [%.3f %.3f; %.3f %.3f] \n", Ad)
fprintf("Qd = [%.3f %.3f; %.3f %.3f] \n", Qd)
fprintf("Rd = %.3f \n", Rd)


%% PART C
fprintf("\n<strong>(c)</strong>\n");

% continuous steady state kalman gain
% P_dot = A*P + P*A' + Bw*Qc*Bw' - P*C'*R^-1*C*P
[Pss, Lss, Sss] = icare(A', C', Bw*Qc*Bw', Rc);

xm = zeros(2,M);
xp = zeros(2,M);
Pm = zeros(2,2,M);
Pp = zeros(2,2,M);
L = zeros(2,M);
% y = zeros(1,M);

% kalman filter
for k = 1:M

    if k == 1
        % initialize
        xm(:,k) = [0; 0];
        Pm(:,:,k) = eye(2);

    else
        % time update (propogate)
        xm(:,k) = Ad * xp(:,k-1);
        Pm(:,:,k) = Ad * Pp(:,:,k-1) * Ad' + Qd;

    end

    % measurement update (correct)
    L(:,k) = Pm(:,:,k) * C' * (C * Pm(:,:,k) * C' + Rd)^-1;
    Pp(:,:,k) = (I - L(:,k) * C) * Pm(:,:,k);
    xp(:,k) = xm(:,k) + L(:,k) * (y(k) - C * xm(:,k));

end

Pp_ss = Pp(:,:,end);
Pm_ss = Pm(:,:,end);
L_ss = L(:,end);

Ad_cl = Ad - L_ss*C;
Poles = eig(Ad_cl);

fprintf("L_ss = [%.3f; %.3f] \n", L_ss)
fprintf("Pp_ss = [%.3f, %.3f; %.3f, %.3f] \n", Pp_ss)
fprintf("Pm_ss = [%.3f, %.3f; %.3f, %.3f] \n", Pm_ss)
fprintf("Ad_cl = [%.3f, %.3f; %.3f, %.3f] \n", Ad_cl)
fprintf("Poles = [%.3f%+.3fi; %.3f%+.3fi] \n", real(Poles(1)), imag(Poles(1)), real(Poles(2)), imag(Poles(2)))

% plot with values
m = 31;
t2 = linspace(0, (m-1)/10, m);

tl = tiledlayout(3,1, TileSpacing="tight", Parent=tab(2));

axes(Parent=tl);
hold on
grid on
plot(t2, L(1,1:m), 'r', LineWidth=2)
plot(t2, L(2,1:m), 'b', LineWidth=2)
legend(sprintf("$L_1$ = %.3f",L_ss(1)), sprintf("$L_2$ = %0.3f",L_ss(2)), Location='eastoutside')
title("\bf{C) Kalman Gain Steady State}")
xlabel("Time [s]")

nexttile;
hold on
grid on
plot(t2, squeeze(Pp(1,1,1:m)), 'r', LineWidth=2)
plot(t2, squeeze(Pp(1,2,1:m)), 'g', LineWidth=2)
plot(t2, squeeze(Pp(2,2,1:m)), 'b', LineWidth=2)
legend(sprintf("$P^+_{1,1}$ = %.3f",Pp_ss(1,1)), sprintf("$P^+_{1,2}$ = %.3f",Pp_ss(1,2)), sprintf("$P^+_{2,2}$ = %.3f",Pp_ss(2,2)), Location='eastoutside')
title("\bf{C) Posteriori Covariance ($P^+$) Steady State}")
xlabel("Time [s]")

nexttile;
hold on
grid on
plot(t2, squeeze(Pm(1,1,1:m)), 'r', LineWidth=2)
plot(t2, squeeze(Pm(1,2,1:m)), 'g', LineWidth=2)
plot(t2, squeeze(Pm(2,2,1:m)), 'b', LineWidth=2)
legend(sprintf("$P^-_{1,1}$ = %.3f",Pm_ss(1,1)), sprintf("$P^-_{1,2}$ = %.3f",Pm_ss(1,2)), sprintf("$P^-_{2,2}$ = %.3f",Pm_ss(2,2)), Location='eastoutside')
title("\bf{C) Priori Covariance ($P^-$) Steady State}")
xlabel("Time [s]")

set(findall(gcf,'-property','FontSize'),'FontSize',16)
set(findall(gcf,'-property','Interpreter'),'Interpreter','latex')

tl = tiledlayout(2,1, Parent=tab(3), TileSpacing="tight");
ax = axes(Parent=tl);
hold on;
grid on;
plot(t, y, '.b', LineWidth=4)
plot(t, xp(1,:), 'g', LineWidth=2)
plot(t, x(1,:), 'r', LineWidth=3)
title("\bf{C) Position Measurements vs Estimates}")
ylabel("Position")

lg = legend("Measurement (y)", "KF Estimate ($\hat{x}^+$)", "State ($x$)", Orientation="horizontal");
lg.Layout.Tile = "north";

nexttile;
hold on;
grid on;
plot(t, xp(2,:), 'g', LineWidth=2)
plot(t, x(2,:), 'r', LineWidth=3)
title("\bf{C) Velocity Estimates}")
ylabel("Position")

xlabel("Time [s]")
set(findall(gcf,'-property','FontSize'),'FontSize',16)
set(findall(gcf,'-property','Interpreter'),'Interpreter','latex')


%% PART D
fprintf("\n<strong>(d)</strong>\n");

% calculate the norm of the standard deviation
%  -> done in PART C :)

N = sqrt(sum(std(x - xp, [], 2).^2, 1));
fprintf("N = %.3f \n", N)


%% PART E
fprintf("\n<strong>(e)</strong>\n");

% changing ratio of Q to R
%  -> keep R constant and change Q

ratio = logspace(-6,6,350);
N = zeros(1,250);

for i = 1:length(ratio)
    R = 1;
    Q = ratio(i) * R;

    % kalman filter
    xm = zeros(2,M);
    xp = zeros(2,M);
    Pm = zeros(2,2,M);
    Pp = zeros(2,2,M);
    L = zeros(2,M);
    for k = 1:M
    
        if k == 1
            % initialize
            xm(:,k) = [0; 0];
            Pm(:,:,k) = eye(2);
    
        else
            % time update (propogate)
            xm(:,k) = Ad * xp(:,k-1);
            Pm(:,:,k) = Ad * Pp(:,:,k-1) * Ad' + Q;
    
        end
    
        % measurement update (correct)
        L(:,k) = Pm(:,:,k) * C' * (C * Pm(:,:,k) * C' + R)^-1;
        Pp(:,:,k) = (I - L(:,k) * C) * Pm(:,:,k);
        xp(:,k) = xm(:,k) + L(:,k) * (y(k) - C * xm(:,k));

    end

    % error standard deviation
    N(i) = sqrt(sum(std(x - xp, [], 2).^2, 1));

end

min_N = min(N);
min_ratio = ratio(N==min_N);

fprintf("Minimim StdDev = %.3f \n", min_N);
fprintf("Ratio = %f \n", min_ratio);

ax = axes(Parent=tab(4));
semilogx(ratio, N, 'b.', LineWidth=4);
hold on;
xline(min_ratio, 'r', LineWidth=1.5, HandleVisibility='off');
grid on;
% text(min_ratio+1e-3, min_N-0.1, 0, ["R:Q = ", num2str(min_ratio), " \sigma = ", num2str(min_N)])
title("\bf{Comparison of Normalized Standard Deviation at Various R:Q Ratios with R=1}")
xlabel("Ratio")
ylabel("$||\sigma||$")
legend(["N (N = "+sprintf("%.3f @ Ratio = %f)",min_N,min_ratio)], Location="northwest")
set(findall(gcf,'-property','FontSize'),'FontSize',16)
set(findall(gcf,'-property','Interpreter'),'Interpreter','latex')

%%

% exportgraphics(tab(1), "./media/p1_a.png");
% exportgraphics(tab(2), "./media/p1_c.png");
% exportgraphics(tab(3), "./media/p1_c2.png");
% exportgraphics(tab(4), "./media/p1_e.png");


