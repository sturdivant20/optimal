%% Optimal HW3 - Problem 4 | Daniel Sturdivant
clc; clear; close all;
fprintf("<strong>PROBLEM 5</strong>\n");

A = [0,0;-1,-2];
B = [1;1];
C = [0,1];
D = 0;

sys = ss(A,B,C,D);
G = tf(sys);

G

figure()
bode(G)

figure()
rlocus(G)