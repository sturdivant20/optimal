%% Daniel Sturdivant || MECH 7710 Optimal || HW2
clc; clear; close all;

% f = figure('units','normalized','position',[0.1 0.1 0.7 0.8]);
f = figure('units','normalized','position',[0.1 0.1 0.9 1.0]);
tabs = uitabgroup(f);
tab(1) = uitab(Parent=tabs, Title='3b)');
tab(2) = uitab(Parent=tabs, Title='3c)');
tab(3) = uitab(Parent=tabs, Title='3d)');
tab(4) = uitab(Parent=tabs, Title='4b)');
tab(5) = uitab(Parent=tabs, Title='4b)-1');
tab(6) = uitab(Parent=tabs, Title='4d)');
tab(7) = uitab(Parent=tabs, Title='4e)');
tab(8) = uitab(Parent=tabs, Title='4e)-1');
tab(9) = uitab(Parent=tabs, Title='5a)');
tab(10) = uitab(Parent=tabs, Title='5b)');

%% PROBLEM 1
fprintf("\n<strong>PROBLEM 1</strong> \n");


%% PROBLEM 2
fprintf("\n<strong>PROBLEM 2</strong> \n");


%% PROBLEM 3
fprintf("\n<strong>PROBLEM 3</strong> \n");

w = rad2deg(2*pi*1); % angular frequency
a = 0.5;    % scale factor
b = 3;      % bias
g = @(t) a*(100*sin(w.*t)) + b + 0.3.*randn(size(t)); % gyro measurement
G = @(t) [100*sin(w.*t), ones(size(t))];              % geometry matrix

% PART A
dt = 0.1; % gyro update rate (1/batch_size)
t = (dt:dt:1)';
y = g(t);
H = G(t);
x = (H'*H)\(H'*y);
fprintf("(a)\nx = [%f %f]'\n", x);

% PART B
fprintf("(b)\n");

dt = 0.1; % gyro update rate / batch size
t0 = 0.1;
x = zeros(2,1000);
P = zeros(2,2,1000);

% 1000 monte carlo runs
for i = 1:1000
    t = (t0:dt:t0+0.9)';
    y = g(t);
    H = G(t);
    x(:,i) = (H'*H)^-1 * (H'*y);
    P(:,:,i) = 0.3^2 * (H'*H)^-1;
    t0 = t0+1;
end

std_est = std(x,[],2);
std_thr = mean(reshape(sqrt([P(1,1,:); P(2,2,:)]), 2,1000),2);
mean_est = mean(x,2);
mean_thr = [0.5;3];
fprintf("sig_est = [%f %f]'\n", std_est);
fprintf("mu_est  = [%f %f]'\n", mean_est);
fprintf("sig_thr = [%f %f]'\n", std_thr);
fprintf("mu_thr  = [%f %f]'\n", mean_thr);

tl = tiledlayout(2,1,Parent=tab(1), TileSpacing="compact");
title(tl(1), "\textbf{Batch Size 10 Monte Carlo}", Interpreter="latex", fontsize=18);
xlabel(tl(1), "Monte Carlo Iteration", fontsize=16);
ylabel(tl(1), "Monte Carlo Estimates", fontsize=16);
ax(1) = axes(Parent=tl);
hold on;
yline(mean_thr(1), Color="#000000", LineWidth=3, DisplayName="Mean");
yline(mean_thr(1)+std_thr(1), Color="#D95319", LineWidth=3, DisplayName="\pm1\sigma");
yline(mean_thr(1)-std_thr(1), Color="#D95319", LineWidth=3, HandleVisibility="off");
plot(1:1000, x(1,:), ".", Color="#4DBEEE", MarkerSize=10, DisplayName="Estimates");
legend;
title("a")
grid on;
nexttile;
hold on;
yline(mean_thr(2), Color="#000000", LineWidth=3, DisplayName="Mean");
yline(mean_thr(2)+std_thr(2), Color="#D95319", LineWidth=3, DisplayName="\pm1\sigma");
yline(mean_thr(2)-std_thr(2), Color="#D95319", LineWidth=3, HandleVisibility="off");
plot(1:1000, x(2,:), ".", Color="#4DBEEE", MarkerSize=10, DisplayName="Estimates");
legend;
title("b")
grid on;
ax(1).XAxis.Visible = "off";

% PART C
fprintf("(c)\n");

dt = 1/1000; % gyro update rate / batch size
t0 = 1/1000;
x = zeros(2,1000);
P = zeros(2,2,1000);

% 1000 monte carlo runs
for i = 1:1000
    t = (t0:dt:t0+999/1000)';
    y = g(t);
    H = G(t);
    x(:,i) = (H'*H)^-1 * (H'*y);
    P(:,:,i) = 0.3^2 * (H'*H)^-1;
    t0 = t0+1;
end

std_est = std(x,[],2);
std_thr = mean(reshape(sqrt([P(1,1,:); P(2,2,:)]), 2,1000),2);
mean_est = mean(x,2);
mean_thr = [0.5;3];
fprintf("sig_est = [%f %f]'\n", std_est);
fprintf("mu_est  = [%f %f]'\n", mean_est);
fprintf("sig_thr = [%f %f]'\n", std_thr);
fprintf("mu_thr  = [%f %f]'\n", mean_thr);

tl = tiledlayout(2,1,Parent=tab(2), TileSpacing="compact");
title(tl, "\textbf{Batch Size 1000 Monte Carlo}", Interpreter="latex", fontsize=18);
xlabel(tl, "Monte Carlo Iteration", fontsize=16);
ylabel(tl, "Monte Carlo Estimates", fontsize=16);
ax(2) = axes(Parent=tl);
hold on;
yline(mean_thr(1), Color="#000000", LineWidth=3, DisplayName="Mean");
yline(mean_thr(1)+std_thr(1), Color="#D95319", LineWidth=3, DisplayName="\pm1\sigma");
yline(mean_thr(1)-std_thr(1), Color="#D95319", LineWidth=3, HandleVisibility="off");
plot(1:1000, x(1,:), ".", Color="#4DBEEE", MarkerSize=10, DisplayName="Estimates");
legend;
title("a")
grid on;
nexttile;
hold on;
yline(mean_thr(2), Color="#000000", LineWidth=3, DisplayName="Mean");
yline(mean_thr(2)+std_thr(2), Color="#D95319", LineWidth=3, DisplayName="\pm1\sigma");
yline(mean_thr(2)-std_thr(2), Color="#D95319", LineWidth=3, HandleVisibility="off");
plot(1:1000, x(2,:), ".", Color="#4DBEEE", MarkerSize=10, DisplayName="Estimates");
legend;
title("b")
grid on;
ax(2).XAxis.Visible = "off";

% PART D
fprintf("(d)\n");

dt = 1/100;
t = (dt:dt:1)';
R = (0.3^2 * eye(1))^-1;
x = zeros(2,100);
P = zeros(2,2,100);

for i = 1:100
    if i == 1
        y = g(t(1));
        H = G(t(1));
        Q = H'*H;
        x(:,1) = Q^-1 * (H'*y);
        P(:,:,1) = 0.3^2 * Q^-1;
    else
        y = g(t(i));
        H = G(t(i));
        Q = Q + H'*H;
        dx = Q^-1 * H' * (y - H*x(:,i-1));

        x(:,i) = x(:,i-1) + dx;
        P(:,:,i) = 0.3^2 * Q^-1;
    end
end

std_thr = reshape(sqrt([P(1,1,:); P(2,2,:)]), 2,100);

tl = tiledlayout(3,1,Parent=tab(3));
title(tl, "\textbf{Recursive Least Squares over 100 Measurements}", Interpreter="latex", fontsize=18);
xlabel(tl, "Time [s]", fontsize=16);
ax(3) = axes(Parent=tl);
hold on;
plot(t, x(1,:), LineWidth=2, DisplayName="$$\hat{a}$$");
plot(t, x(2,:), LineWidth=2, DisplayName="$$\hat{b}$$");
ylabel("Estimates")
legend(Interpreter="latex");
grid on;
ax(4) = nexttile;
hold on;
plot(t, std_thr(1,:), LineWidth=2, DisplayName="$$\sigma_a$$");
plot(t, std_thr(2,:), LineWidth=2, DisplayName="$$\sigma_b$$");
ylabel("Standard Deviation");
legend(Interpreter="latex");
ylim([0,1]);
grid on;
nexttile;
hold on;
plot(t, x(1,:)-0.5, LineWidth=2, DisplayName="Error in $$\hat{a}$$");
plot(t, x(2,:)-3, LineWidth =2, DisplayName="Error in $$\hat{b}$$");
ylabel("Error");
legend(Interpreter="latex");
grid on;
ax(3).XAxis.Visible = "off";
ax(4).XAxis.Visible = "off";
set(findall(gcf, '-property','FontSize'), FontSize=14);


%% PROBLEM 4
fprintf("\n<strong>PROBLEM 4</strong> \n");

% PART B
fprintf("(b)\n");

% discrete system
num = 0.25 .* [1, -0.8]; 
den = [1, -1.9, 0.95];
G = tf(num,den, 1);

% simulated discrete system
U = randn(1000,1);
y = dlsim(num,den,U); 
sigma = 0.01;
Y = y + sigma*randn(1000,1); 

% least squares coefficient estimation
H = [U(2:end-1), U(1:end-2), -Y(2:end-1), -Y(1:end-2)];
m = Y(3:end);
x = (H'*H)^-1*H'*m;
P = sigma^2 * (H'*H)^-1;

% recreated system
numid = [x(1), x(2)];
denid = [1, x(3), x(4)];
Gid = tf(numid,denid, 1);
SNR = std(Y) / sigma;

fprintf("x   = [%f, %f, %f, %f]'\n", x);
fprintf("SNR = %f\n", SNR);

tl = tiledlayout(2,1,Parent=tab(4), TileSpacing="tight");
title(tl, "\textbf{System ID Bode Plot}", Interpreter="latex", fontsize=18);
ax(5) = axes(Parent=tl);
[m, p, w] = bode(G);
[mid, pid, wid] = bode(Gid);
semilogx(w,20*log10(m(:)), LineWidth=2, DisplayName="Original");
hold on;
semilogx(wid,20*log10(mid(:)), LineWidth=2, DisplayName="Estimated");
xline(w(end)+0.01, Color="k", LineWidth=2, HandleVisibility="off");
grid on;
legend;
ylabel("Magnitude [dB]");
ax(6) = nexttile;
semilogx(w,p(:), LineWidth=2, DisplayName="Original");
hold on;
semilogx(wid,pid(:), LineWidth=2, DisplayName="Estimated");
xline(w(end)+0.01, Color="k", LineWidth=2, HandleVisibility="off");
grid on;
legend;
ylabel("Phase [deg]");
xlabel("Frequency [rad/s]");
ax(5).XAxis.Visible = "off";
ax(5).FontSize = 14;
ax(6).FontSize = 14;
linkaxes([ax(5), ax(6)], 'x');

tl = tiledlayout(1,1,Parent=tab(5), TileSpacing="tight");
title(tl, "\textbf{Signal Plot}", Interpreter="latex", fontsize=18);
ax(7) = axes(Parent=tl);
hold on;
plot(y, LineWidth=2, DisplayName="y");
plot(Y, LineWidth=2, DisplayName="Y");
grid on;
legend;
xlabel("Samples");
ylabel("Magnitude");
ax(7).FontSize = 14;

% PART C
fprintf("(c)\n");

U = zeros(1000,10);
Y = zeros(1000,10);
x = zeros(4,10);
P = zeros(4,4,10);
for i = 1:10
    % resimulate system
    U(:,i) = randn(1000,1);
    y = dlsim(num,den,U(:,i)); 
    sigma = 0.01;
    Y(:,i) = y + sigma*randn(1000,1); 
    
    % least squares
    H = [U(2:end-1,i), U(1:end-2,i), -Y(2:end-1,i), -Y(1:end-2,i)];
    m = Y(3:end,i);
    x(:,i) = (H'*H)^-1*H'*m;
    P(:,:,i) = sigma^2 * (H'*H)^-1;
end

std_est = std(x,[],2);  % theoretical = [0.01,  0.01,  0.01, 0.01]
mean_est = mean(x,2);   % theoretical = [0.25, -0.20, -1.90, 0.95]
std_thr = mean(reshape(sqrt([P(1,1,:); P(2,2,:); P(3,3,:); P(4,4,:)]), 4,10),2);
mean_thr = [0.25, -0.20, -1.90, 0.95];

fprintf("sig_est = [%f, %f, %f, %f]'\n", std_est);
fprintf("mu_est  = [%f, %f, %f, %f]'\n", mean_est);
fprintf("sig_thr = [%f, %f, %f, %f]'\n", std_thr);
fprintf("mu_thr  = [%f, %f, %f, %f]'\n", mean_thr);

% PART D
sigma = 0.1 : 0.1 : 1;
L = length(sigma);
LL = 250;
U = zeros(1000,LL,L);
Y = zeros(1000,LL,L);
x = zeros(4,LL,L);
P = zeros(4,4,LL,L);
std_est = zeros(4,L);
mean_est = zeros(4,L);
SNR = zeros(1,L);

% 10 monte carlo runs
for s = 1:L
    % standard deviations between 0.1 and 1
    for i = 1:LL
        % resimulate system
        U(:,i,s) = randn(1000,1);
        y = dlsim(num,den,U(:,i)); 
        Y(:,i,s) = y + sigma(s)*randn(1000,1); 
        
        % least squares
        H = [U(2:end-1,i,s), U(1:end-2,i,s), -Y(2:end-1,i,s), -Y(1:end-2,i,s)];
        m = Y(3:end,i,s);
        x(:,i,s) = (H'*H)^-1*H'*m;
        P(:,:,i,s) = sigma(s)^2 * (H'*H)^-1;
    end
    std_est(:,s) = std(x(:,:,s), [], 2);
    mean_est(:,s) = mean(x(:,:,s), 2);
    SNR(:,s) = std(Y(:,:,s),[],'all') / sigma(s);
end

tl = tiledlayout(3,1,Parent=tab(6), TileSpacing="tight");
title(tl, "\textbf{Varying $\sigma$}", Interpreter="latex", fontsize=18);
ax(8) = axes(Parent=tl);
hold on;
% plot(sigma, sigma, 'k--', HandleVisibility="off");
plot(sigma, std_est(1,:), Color="#0072BD", LineWidth=2, DisplayName="\sigma_A");
plot(sigma, std_est(2,:), Color="#77AC30", LineWidth=2, DisplayName="\sigma_B");
plot(sigma, std_est(3,:), Color="#EDB120", LineWidth=2, DisplayName="\sigma_C");
plot(sigma, std_est(4,:), Color="#A2142F", LineWidth=2, DisplayName="\sigma_D");
grid on;
legend;
ylabel("\sigma_{Y}");
ax(9) = nexttile;
hold on;
plot(sigma, mean_est(1,:), Color="#0072BD", LineWidth=2, DisplayName="\mu_A");
plot(sigma, mean_est(2,:), Color="#77AC30", LineWidth=2, DisplayName="\mu_B");
plot(sigma, mean_est(3,:), Color="#EDB120", LineWidth=2, DisplayName="\mu_C");
plot(sigma, mean_est(4,:), Color="#A2142F", LineWidth=2, DisplayName="\mu_D");
grid on;
legend;
ylabel("\mu_Y");
ax(10) = nexttile;
plot(sigma, SNR, LineWidth=2, DisplayName="SNR");
grid on;
legend;
ylabel("SNR")
xlabel("\sigma");
ax(8).XAxis.Visible = "off";
ax(9).XAxis.Visible = "off";
ax(8).FontSize = 14;
ax(9).FontSize = 14;
ax(10).FontSize = 14;
% linkaxes([ax(8), ax(9), ax(10)], 'x');


Gid = tf([mean_est(1,10), mean_est(2,10)], [1, mean_est(3,10), mean_est(4,10)], 1);
tl = tiledlayout(2,1,Parent=tab(7), TileSpacing="tight");
title(tl, "\textbf{System ID Bode Plot}", Interpreter="latex", fontsize=18);
ax(11) = axes(Parent=tl);
[m, p, w] = bode(G);
[mid, pid, wid] = bode(Gid);
semilogx(w,20*log10(m(:)), LineWidth=2, DisplayName="Original");
hold on;
semilogx(wid,20*log10(mid(:)), LineWidth=2, DisplayName="Estimated");
xline(w(end)+0.01, Color="k", LineWidth=2, HandleVisibility="off");
grid on;
legend;
ylabel("Magnitude [dB]");
ax(12) = nexttile;
semilogx(w,p(:), LineWidth=2, DisplayName="Original");
hold on;
semilogx(wid,pid(:), LineWidth=2, DisplayName="Estimated");
xline(w(end)+0.01, Color="k", LineWidth=2, HandleVisibility="off");
grid on;
legend;
ylabel("Phase [deg]");
xlabel("Frequency [rad/s]");
ax(11).XAxis.Visible = "off";
ax(11).FontSize = 14;
ax(12).FontSize = 14;
% linkaxes([ax(11), ax(12)], 'x');

tl = tiledlayout(1,1,Parent=tab(8), TileSpacing="tight");
title(tl, "\textbf{Signal Plot}", Interpreter="latex", fontsize=18);
ax(13) = axes(Parent=tl);
hold on;
plot(y, LineWidth=2, DisplayName="y");
plot(mean(Y(:,:,10),2), LineWidth=2, DisplayName="Y");
grid on;
legend;
xlabel("Samples");
ylabel("Magnitude");
ax(13).FontSize = 14;


%% PROBLEM 5
fprintf("\n<strong>PROBLEM 5</strong> \n");
Tw = 1;
wc = pi;
A = 1;
Sy = @(x) A ./ (Tw^2.*x.^2 + 1);


G = tf([1], [Tw, 1]);
[m_out, p_out, w_out] = bode(G);

tl = tiledlayout(3,1, TileSpacing="tight", Parent=tab(9));
xlabel(tl, "Angular Frequency", FontSize=16);
ylabel(tl, "Magnitude", FontSize=16);
ax1 = axes(Parent=tl);
hold on;
yline(0, 'k', LineWidth=1.5);
xline(0, 'k', LineWidth=1.5);
xline(wc, '--k', LineWidth=1.5);
xline(-wc, '--k', LineWidth=1.5);
plot([-wc, wc], [A, A], 'r', LineWidth=2.5);
grid on;
xlim([-10,10]);
ylim([-1,2]);
xticks([-wc, wc]);
xticklabels({'-\omega_c', '\omega_c'});
yticks([A]);
yticklabels({'A'});
title("Case 1");
ax2 = nexttile;
hold on;
yline(0, 'k', LineWidth=1.5);
xline(0, 'k', LineWidth=1.5);
xline(wc, '--k', LineWidth=1.5);
xline(-wc, '--k', LineWidth=1.5);
plot([-10, 10], [A, A], 'r', LineWidth=2.5);
grid on;
xlim([-10,10]);
ylim([-1,2]);
xticks([-wc, wc]);
xticklabels({'-\omega_c', '\omega_c'});
yticks([A]);
yticklabels({'A'});
title("Case 2");
ax3 = nexttile;
semilogx(w_out, 20*log10(m_out(:))+Tw, 'r', LineWidth=2.5);
hold on;
yline(0, 'k', LineWidth=1.5);
xline(0, 'k', LineWidth=1.5);
xline(1/Tw, '--k', LineWidth=1.5);
grid on;
xticks([1/Tw]);
xticklabels({'1/T_{\omega}'});
yticks([Tw]);
yticklabels({'T_{\omega}'});
% ylim([-1,2]);
title("|G(j\omega)|");
set(findall(gcf, '-property','FontSize'), FontSize=16);

tl = tiledlayout(2,1, TileSpacing="tight", Parent=tab(10));
xlabel(tl, "Angular Frequency", FontSize=16);
ylabel(tl, "Magnitude", FontSize=16);
ax1 = axes(Parent=tl);
x1 = linspace(-wc/8.5,wc/8.5, 100);
hold on;
plot(x1, Sy(x1), 'r', LineWidth=2.5);
symlog(gca, 'xy', -0.95)
grid off;
yline(0, 'k', LineWidth=1.5);
xline(0, 'k', LineWidth=1.5);
xline(wc/5, '--k', LineWidth=1.5);
xline(-wc/5, '--k', LineWidth=1.5);
xline(1/Tw/5, '--k', LineWidth=1.5);
xline(-1/Tw/5, '--k', LineWidth=1.5);
xlim([-1.95,1.95]);
% ylim([-1,2]);
xticks([-wc/5, -1/Tw/5, 1/Tw/5 ,wc/5]);
yticks([A]);
xticklabels({'-\omega_c', '-1/Tw', '1/Tw', '\omega_c'});
yticklabels({'A'});
ax1.XAxis.TickLength = [0 0];
ax1.YAxis.TickLength = [0 0];
ax2 = nexttile;
x2 = linspace(-10,10, 100);
hold on;
plot(x2, Sy(x2), 'r', LineWidth=2.5);
symlog(gca, 'xy', -0.95)
grid off;
yline(0, 'k', LineWidth=1.5);
xline(0, 'k', LineWidth=1.5);
xline(wc/5, '--k', LineWidth=1.5);
xline(-wc/5, '--k', LineWidth=1.5);
xline(1/Tw/5, '--k', LineWidth=1.5);
xline(-1/Tw/5, '--k', LineWidth=1.5);
xlim([-1.95,1.95]);
% ylim([-1,2]);
xticks([-wc/5, -1/Tw/5, 1/Tw/5, wc/5]);
yticks([A]);
xticklabels({'-\omega_c', '-1/Tw', '1/Tw', '\omega_c'});
yticklabels({'A'});
ax2.XAxis.TickLength = [0 0];
ax2.YAxis.TickLength = [0 0];
set(findall(gcf, '-property','FontSize'), FontSize=16);


%%
% exportgraphics(tab(1), "./media/3b.png");
% exportgraphics(tab(2), "./media/3c.png");
% exportgraphics(tab(3), "./media/3d.png");
% exportgraphics(tab(4), "./media/4b.png");
% exportgraphics(tab(5), "./media/4b-1.png");
% exportgraphics(tab(6), "./media/4d.png");
% exportgraphics(tab(7), "./media/4e.png");
% exportgraphics(tab(8), "./media/4e-1.png");
exportgraphics(tab(9), "./media/5a.png");
exportgraphics(tab(10), "./media/5b.png");
