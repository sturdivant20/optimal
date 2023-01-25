%% Daniel Sturdivant || MECH 7710 - Optimal || HW0
clc; close all; clear;

f = figure(Name="Optimal HW0");
tabs = uitabgroup(f);
tab1 = uitab(tabs, Title="Observer");
tab2 = uitab(tabs, Title="Controller + Observer");
tab3 = uitab(tabs, Title="Comp. Bode Plot");
tab4 = uitab(tabs, Title="Discrete Poles");
tab5 = uitab(tabs, Title="Cont. vs Disc.");

%% Problem 1
fprintf("PROBLEM 1 \n");

syms z1 z2 z1d z2d u y
J = 10; % [kg-m^2]
b = 1;  % [N-m-s/rad]

f = [z2; 1/J*(b*-z2 + u)]; 
g = [z1];

A = double(jacobian(f, [z1, z2]));
B = double(jacobian(f, [u]));
C = double(jacobian(g, [z1, z2]));
D = double(jacobian(g, [u]));

s = eig(A);
fprintf("The open-loop eigenvales are s = %.2f, %.2f \n", s);
sys = ss(A,B,C,D);
TF = tf(sys);

fprintf("\n");

%% Problem 2
fprintf("PROBLEM 2 \n");

Ob = obsv(A,C);
if rank(Ob) == size(A,1)
    fprintf("System is observable, rank(Ob) = %.0f \n", rank(Ob));
else
    fprintf("System is unobservable, rank(Ob) = %.0f \n", rank(Ob));
end

f = 50; % [Hz]
l = 0.7;
w = 2*pi*f; % [rad/s]

syms s
sL = double(solve(s^2 + 2*w*l*s + w^2, s));
fprintf("The poles when f = %.2f and l = %.2f are s = %.2f%+.2fi, %.2f%+.2fi \n", f, l, real(sL(1)), imag(sL(1)), real(sL(2)), imag(sL(2)));

L = place(A', C', sL)';
fprintf("L = [%.2f \n     %.2f] \n", L)

% closed-loop system
sys_L = ss(A-L*C, B, C, D);

axes(Parent=tab1);
[y_obsv,t] = step(sys_L, 0.2);
plot(t, squeeze(y_obsv), 'b', LineWidth=3);
title("2) Step Response of Observer", FontSize=18);
xlabel('time [s]', FontSize=16);
ylabel('Amplitude', FontSize=16);

fprintf("\n");

%% PROBLEM 3
fprintf("PROBLEM 3 \n");

Co = ctrb(A,B);
if rank(Co) == size(A,1)
    fprintf("System is controllable, rank(Co) = %.0f \n", rank(Co));
else
    fprintf("System is uncontrollable, rank(Co) = %.0f \n", rank(Co));
end

f = 10; % [Hz]
l = 0.7;
w = 2*pi*f; % [rad/s]

syms s
sK = double(solve(s^2 + 2*w*l*s + w^2, s));
fprintf("The poles when f = %.2f and l = %.2f are s = %.2f%+.2fi, %.2f%+.2fi \n", f, l, real(sK(1)), imag(sK(1)), real(sK(2)), imag(sK(2)));

K = place(A, B, sK);
fprintf("K = [%.2f, %.2f] \n", K)

% closed loop system (controller + observer)
Acl = [A-B*K         , B*K; ...
       zeros(size(A)), A-L*C];
Bcl = [B; ...
       zeros(size(B))];
Ccl = [C, zeros(size(C))];
Dcl = D;
sys_K = ss(Acl, Bcl, Ccl, Dcl);

axes(Parent=tab2);
[y_cont, t] = step(sys_K, 0.2);
plot(t, y_cont, 'r', LineWidth=3);
title("2) Step Response of Controller", FontSize=18);
xlabel('time [s]', FontSize=16);
ylabel('Amplitude', FontSize=16);

fprintf("\n");

%% PROBLEM 4
fprintf("PROBLEM 4 \n");

A_comp = A - B*K - L*C;
B_comp = L;
C_comp = -K;

sys_comp = ss(A_comp, B_comp, C_comp, 0);
TF_comp = tf(sys_comp)

% f_path = TF_comp * TF;
% TF_CL = feedback(TF*TF_comp,1);
syms s
TF_CL = -C*(s*eye(size(A_comp)) - A_comp)^-1*B;
[n,d] = numden(simplify(TF_CL));
d = flip(double(coeffs(d)));
n = flip(double(coeffs(n))) / d(1);
d = d / d(1);
TF_CL = tf(n, d)

[Gm,Pm] = margin(TF);
fprintf("Open-Loop: Gain Margin: %f, Phase Margin: %f \n", 20*log10(Gm), Pm);
[Gm,Pm] = margin(TF_CL);
fprintf("Closed-Loop: Gain Margin: %f, Phase Margin: %f \n", 20*log10(Gm), Pm);

axes(Parent=tab3);
sgtitle('3) Bode Plot', FontSize=18, FontWeight='Bold')
w = logspace(-2, 6, 10000);
[mag_CL, phase_CL, w_CL] = bode(TF_CL, w);
subplot(2,1,1);
semilogx(w_CL, 20*log10(squeeze(mag_CL)), 'b', LineWidth=3);
yline(0, LineWidth=2);
ylabel('Magnitude [dB]', FontSize=16);
grid on;
subplot(2,1,2);
semilogx(w_CL, squeeze(phase_CL), 'r', LineWidth=3);
yline(0, LineWidth=2);
xlabel('Frequency [rad/s]', FontSize=16);
ylabel('Phase [deg]', FontSize=16);
grid on;


% axes(Parent=tab3);
% % figure();
% sgtitle('3) Bode Plot', FontSize=18, FontWeight='Bold')
% w = logspace(-4, 5, 10000);
% [mag_OL, phase_OL, w_OL] = bode(TF, w);
% [mag_CL, phase_CL, w_CL] = bode(TF_CL, w);
% [mag_, phase_, w_] = bode(f_path, w);
% mag_OL = 20*log10(squeeze(mag_OL));
% mag_CL = 20*log10(squeeze(mag_CL));
% mag_ = 20*log10(squeeze(mag_));
% phase_OL = squeeze(phase_OL);
% phase_CL = squeeze(phase_CL);
% phase_ = squeeze(phase_);
% 
% subplot(2,1,1);
% semilogx(w_OL, mag_OL, 'b', LineWidth=3);
% hold on;
% semilogx(w_, mag_, 'g', LineWidth=3);
% semilogx(w_CL, mag_CL, 'r', LineWidth=3);
% yline(0);
% hold off;
% xlim([1e-4,1e5]);
% ylabel('Magnitude [db]', FontSize=16);
% grid on;
% legend({'Uncompensated', 'Forward Path', 'Compensated'}, FontSize=16)
% subplot(2,1,2);
% semilogx(w_OL, phase_OL, 'b', LineWidth=3);
% hold on;
% semilogx(w_, phase_, 'g', LineWidth=3);
% semilogx(w_CL, phase_CL, 'r', LineWidth=3);
% yline(0);
% hold off;
% xlim([1e-4,1e5]);
% xlabel('Frequency [rad/s]', FontSize=16);
% ylabel('Phase [deg]', FontSize=16);
% grid on;
% % 
% [Gm,Pm] = margin(TF);
% fprintf("Open-Loop: Gain Margin: %f, Phase Margin: %f \n", Gm, Pm);
% [Gm,Pm] = margin(TF_CL);
% fprintf("Closed-Loop: Gain Margin: %f, Phase Margin: %f \n", Gm, Pm);
% [Gm, Pm] = margin(f_path);
% fprintf("Foward-Path: Gain Margin: %f, Phase Margin: %f \n", Gm, Pm);

fprintf("\n")

%% PROBLEM 5
fprintf("PROBLEM 5 \n");

sys2 = c2d(sys, 1/1000);
TF2 = tf(sys2);
Az = sys2.A;
Bz = sys2.B;
Cz = sys2.C;
Dz = sys2.D;
z = eig(Az);
fprintf("The open-loop eigenvales are z = %.4f, %.4f \n", z);

zL = exp(sL*1/1000);
zK = exp(sK*1/1000);
Lz = place(Az', Cz', zL)';
Kz = place(Az, Bz, zK);
fprintf("L = [%.2f \n     %.2f] \n", L)
fprintf("K = [%.2f, %.2f] \n", K)

A_c2 = Az - Bz*Kz - Lz*Cz;
B_c2 = Lz;
C_c2 = -Kz;

syms z
TF_comp2 = -Kz*(z*eye(size(A_c2)) - A_c2)^-1*Lz;
[n,d] = numden(simplify(TF_comp2));
d = flip(double(coeffs(d)));
n = flip(double(coeffs(n))) / d(1);
d = d / d(1);
TF_comp2 = tf(n, d, 1/1000)

syms z
TF_CL2 = -Cz*(z*eye(size(A_c2)) - A_c2)^-1*Bz;
[n,d] = numden(simplify(TF_CL2));
d = flip(double(coeffs(d)));
n = flip(double(coeffs(n))) / d(1);
d = d / d(1);
TF_CL2 = tf(n, d, 1/1000)

axes(Parent=tab4);
title('5) Discrete Poles', FontSize=18);
hold on;
plot(real(zL), imag(zL), 'bx', MarkerSize=10);
plot(real(zK), imag(zK), 'rx', MarkerSize=10);
% unit circle
ang = 0:0.01:2*pi; xp = cos(ang); yp = sin(ang);
plot(xp, yp, 'k', LineWidth=3);
yline(0, LineWidth=2);
xline(0, LineWidth=2);
hold off;
grid on;
legend({'Observer Poles', 'Controller Poles'}, FontSize=16);
xlabel('Real', FontSize=16);
ylabel('Imaginary', FontSize=16);

fprintf("\n");

%% PROBLEM 5
fprintf("PROBLEM 5 \n");

[Ys,t1] = step(TF_comp, 0.2);
[Yz,t2] = step(TF_comp2, 0.2);

axes(Parent=tab5)
title('Continuous vs. Discrete Compensator', FontSize=18)
hold on;
plot(t1, Ys, 'b', LineWidth=3);
plot(t2, Yz, 'r', LineWidth=3);
hold off;
grid on;
legend({'Continuous', 'Discrete'}, FontSize=16);
xlabel('Time [s]', FontSize=16);
ylabel('Amplitude', FontSize=16);

fprintf("\n");
