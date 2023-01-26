%% Daniel Sturdivant || MECH 7710 - Optimal || HW0
% https://web.mit.edu/16.31/www/Fall06/1631_topic17.pdf
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

% open loop system
f = [z2; 1/J*(b*-z2 + u)]; 
g = [z1];

A = double(jacobian(f, [z1, z2]));
B = double(jacobian(f, [u]));
C = double(jacobian(g, [z1, z2]));
D = double(jacobian(g, [u]));

% eignevalues and characteristic equation
s = eig(A);
fprintf("The open-loop eigenvales are s = %.2f, %.2f \n", s);
sys = ss(A,B,C,D);
TF = tf(sys);

fprintf("\n");

%% Problem 2
fprintf("PROBLEM 2 \n");

% observability
Ob = obsv(A,C);
if rank(Ob) == size(A,1)
    fprintf("System is observable, rank(Ob) = %.0f \n", rank(Ob));
else
    fprintf("System is unobservable, rank(Ob) = %.0f \n", rank(Ob));
end

% desired observer dynamics
f = 50; % [Hz]
l = 0.7;
w = 2*pi*f; % [rad/s]

syms s
sL = double(solve(s^2 + 2*w*l*s + w^2, s));
fprintf("The poles when f = %.2f and l = %.2f are s = %.2f%+.2fi, " + ...
    "%.2f%+.2fi \n", f, l, real(sL(1)), imag(sL(1)), real(sL(2)), ...
    imag(sL(2)));

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

% controllability
Co = ctrb(A,B);
if rank(Co) == size(A,1)
    fprintf("System is controllable, rank(Co) = %.0f \n", rank(Co));
else
    fprintf("System is uncontrollable, rank(Co) = %.0f \n", rank(Co));
end

% desired  controller dynamics
f = 10; % [Hz]
l = 0.7;
w = 2*pi*f; % [rad/s]

syms s
sK = double(solve(s^2 + 2*w*l*s + w^2, s));
fprintf("The poles when f = %.2f and l = %.2f are s = %.2f%+.2fi, " + ...
    "%.2f%+.2fi \n", f, l, real(sK(1)), imag(sK(1)), real(sK(2)), ...
    imag(sK(2)));

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

% compensator
A_comp = A - B*K - L*C;
B_comp = L;
C_comp = K;

syms s
TF_comp = -K*(s*eye(size(A_comp)) - A_comp)^-1*L;
[n,d] = numden(simplify(TF_comp));
d = flip(double(coeffs(d)));
n = flip(double(coeffs(n))) / d(1);
d = d / d(1);
TF_comp = tf(n, d)

% open loop dynamics
% TF_loop_dynamics = TF_comp * TF;

% closed loop dynamics
A_comp_cl = [A       , B*C_comp; ...
            -B_comp*C, A_comp];
B_comp_cl = [zeros(size(B)); ... 
             B_comp];
C_comp_cl = [C, zeros(size(C_comp))];
sys_comp_CL = ss(A_comp_cl, B_comp_cl, C_comp_cl, 0);
TF_comp_CL = tf(sys_comp_CL)


% gain and phase margins
% [Gm,Pm] = margin(TF);
% fprintf("Open-Loop: Gain Margin: %f, Phase Margin: %f \n", 20*log10(Gm), Pm);
[Gm,Pm] = margin(TF_comp_CL);
fprintf("Closed-Loop: Gain Margin: %f, Phase Margin: %f \n", 20*log10(Gm), Pm);

% bode plot
w = logspace(-2, 6, 10000);
[mag_CL, phase_CL, w_CL] = bode(TF_comp_CL, w);

% [mag_OL, phase_OL, w_OL] = bode(TF_comp, w);
% [mag, phase, w] = bode(sys, w);

axes(Parent=tab3);
sgtitle('3) Bode Plot', FontSize=18, FontWeight='Bold')
subplot(2,1,1);
semilogx(w_CL, 20*log10(squeeze(mag_CL)), 'r', LineWidth=3);
hold on
% semilogx(w_OL, 20*log10(squeeze(mag_OL)), 'b', LineWidth=3);
% semilogx(w, 20*log10(squeeze(mag)), 'g', LineWidth=3);
% legend({'Comp CL', 'Comp', 'OL'})
yline(0, LineWidth=2);
ylabel('Magnitude [dB]', FontSize=16);
grid on;
subplot(2,1,2);
semilogx(w_CL, squeeze(phase_CL), 'r', LineWidth=3);
hold on
% semilogx(w_OL, squeeze(phase_OL), 'b', LineWidth=3);
% semilogx(w, squeeze(phase), 'g', LineWidth=3);
% legend({'Comp CL', 'Comp', 'OL'})
yline(0, LineWidth=2);
xlabel('Frequency [rad/s]', FontSize=16);
ylabel('Phase [deg]', FontSize=16);
grid on;

fprintf("\n")

%% PROBLEM 5
fprintf("PROBLEM 5 \n");

% discrete system
sys2 = c2d(sys, 1/1000);
TF2 = tf(sys2);
Az = sys2.A;
Bz = sys2.B;
Cz = sys2.C;
Dz = sys2.D;
z = eig(Az);
fprintf("The open-loop eigenvales are z = %.4f, %.4f \n", z);

% discrete pole placement
zL = exp(sL*1/1000);
zK = exp(sK*1/1000);
fprintf("Discrete observer poles = %.2f%+.2fi, " + ...
    "%.2f%+.2fi \n",real(zL(1)), imag(zL(1)), real(zL(2)), ...
    imag(zL(2)));
fprintf("Discrete controller poles = %.2f%+.2fi, " + ...
    "%.2f%+.2fi \n",real(zK(1)), imag(zK(1)), real(zK(2)), ...
    imag(zK(2)));
Lz = place(Az', Cz', zL)';
Kz = place(Az, Bz, zK);
fprintf("L = [%.4f \n     %.4f] \n", Lz)
fprintf("K = [%.4f, %.4f] \n", Kz)

% discrete compensator
Az_comp = Az - Bz*Kz - Lz*Cz;
Bz_comp = Lz;
Cz_comp = Kz;

syms z
TF_comp2 = -Kz*(z*eye(size(Az_comp)) - Az_comp)^-1*Lz;
[n,d] = numden(simplify(TF_comp2));
d = flip(double(coeffs(d)));
n = flip(double(coeffs(n))) / d(1);
d = d / d(1);
TF_comp2 = tf(n, d, 1/1000)

% closed loop dynamics
% A_comp_cl2 = [Az         , Bz*Cz_comp; ...
%               -Bz_comp*Cz, Az_comp];
% B_comp_cl2 = [zeros(size(B)); ... 
%              Bz_comp];
% C_comp_cl2 = [Cz, zeros(size(Cz_comp))];
% sys_comp_CL2 = ss(A_comp_cl2, B_comp_cl2, C_comp_cl2, 0);
% TF_comp_CL2 = tf(sys_comp_CL2)
sys_comp_CL2 = c2d(sys_comp_CL, 1/1000);
TF_comp_CL2 = c2d(TF_comp_CL, 1/1000)

axes(Parent=tab4);
title('5) Discrete Poles', FontSize=18);
hold on;
plot(real(zL), imag(zL), 'bx', MarkerSize=10, LineWidth=3);
plot(real(zK), imag(zK), 'rx', MarkerSize=10, LineWidth=3);
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

%% PROBLEM 6
fprintf("PROBLEM 6 \n");

% figure
% step(sys_comp_CL, 0.2, 'b', LineWidth=3);
% hold on
% step(sys_comp_CLz, 0.2, 'r', LineWidth=3);
% hold off

[Ys,t1] = step(sys_comp_CL, 0:0.001:0.2);
[Yz,t2] = step(sys_comp_CL2, 0:0.001:0.2);

axes(Parent=tab5)
title('Continuous vs. Discrete Compensator', FontSize=18)
hold on;
plot(t1, Ys, 'b', LineWidth=3);
stairs(t2, squeeze(Yz), 'r', LineWidth=3);
hold off;
grid on;
legend({'Continuous', 'Discrete'}, FontSize=16);
xlabel('Time [s]', FontSize=16);
ylabel('Amplitude', FontSize=16);

fprintf("\n");
