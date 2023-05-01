%% Optimal HW4 - Problem 2 | Daniel Sturdivant
clc; close all; clear;
fprintf("<strong>PROBLEM 2</strong>\n");

f = figure;
tbs = uitabgroup(Parent=f);
tab(1) = uitab(Parent=tbs, Title="SysID");
tab(2) = uitab(Parent=tbs, Title="A");
tab(3) = uitab(Parent=tbs, Title="B");

G = tf(0.25.*[1, -0.8], [1, -1.9, 0.95], 1);
[m,g,p] = bode(G);

% resimulate system
U = randn(1000,1);
y = dlsim(0.25.*[1, -0.8], [1, -1.9, 0.95], U); 
sigma = 0.01;
Y = y + sigma*randn(1000,1);

sys = arx([Y U], [2 2 1]);
[m_,g_,p_] = bode(sys);

tl = tiledlayout(2,1, Parent=tab(1), TileSpacing="tight");
title(tl, "Bode Plot")
axes(Parent=tl);
semilogx(p, 20*log10(m(:)));
hold on;
semilogx(p_, 20*log10(m_(:)));
grid on;
ylabel("Gain [dB]");
nexttile;
semilogx(p, g(:));
hold on;
semilogx(p_, g_(:));
grid on;
ylabel("Phase [deg]");
xlabel("Frequency [rad/s]");


%% PART A
fprintf("<strong>(a)</strong>\n");

% resimulate system
U = randn(1000,1);
y = dlsim(0.25.*[1, -0.8], [1, -1.9, 0.95], U); 
sigma = 1;
Y = y + sigma*randn(1000,1);

sys_a = arx([Y U], [11 11 1]);
[ma,ga,pa] = bode(sys_a);

tl = tiledlayout(2,1, Parent=tab(2), TileSpacing="tight");
title(tl, "Bode Plot")
axes(Parent=tl);
semilogx(p, 20*log10(m(:)));
hold on;
semilogx(pa, 20*log10(ma(:)));
grid on;
ylabel("Gain [dB]");
nexttile;
semilogx(p, g(:));
hold on;
semilogx(pa, ga(:));
grid on;
ylabel("Phase [deg]");
xlabel("Frequency [rad/s]");


%% PART B
fprintf("<strong>(b)</strong>\n");

