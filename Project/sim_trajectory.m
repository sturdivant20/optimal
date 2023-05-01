clc; clear; close all;

%% TRAJECTORY

Ns = 2;
dt = 1;
toa = [0, 10, 22, 30, 35, 40, 45, 50, 58, 70];

% theta = 0 : 2*pi/8 : 2*pi;
% x = cos(theta);
% y = sin(theta);
% wp = [1,x;
%       0,y] .* 60;
wp = repmat([1;0], 1,length(toa)) .* 80;
wp2 = -wp;

[pos,vel,acc,jerk,pp,~,t] = minjerkpolytraj(wp, toa, 71);
[pos2,vel2,acc2,jerk2,pp2,~,t2] = minjerkpolytraj(wp2, toa, 71);
pos = [pos; pos2];
vel = [vel; vel2];

%% OBSERVATIONS

x_sensor = [-100, -100; ...
               0, -100; 
             100, -100];
% x_sensor = [-100, -500; ...
%                0, -500; 
%              100, -500];
% x_sensor = [-90, -60; ...
%              60,  90; 
%              90, -90;
%             -60,  90].*1;

j = 1;
y = zeros(Ns*size(x_sensor,1), 71);
for i = 1:2:size(pos,1)
    y(j:j+2,:) = wrapTo360(atan2d(x_sensor(:,2) - pos(i+1,:), x_sensor(:,1) - pos(i,:)));
    j = j+3;
end

%% PLOTTING

figure;
hold on;
plot(x_sensor(:,1), x_sensor(:,2), 'bo', LineWidth=3);
for i = 1:2:size(pos,1)
    plot(pos(i,:), pos(i+1,:));
    plot(pos(i,1), pos(i+1,1), 'rx', LineWidth=3);
end
j = 1;
for i = 1:size(y,1)
    j = mod(i,size(x_sensor,1));
    if j == 0
        j = size(x_sensor,1);
    end
    plot([x_sensor(j,1), x_sensor(j,1) - 1000*cosd(y(i,1))], ...
         [x_sensor(j,2), x_sensor(j,2) - 1000*sind(y(i,1))], 'k');
end
xlim([-500, 500]);
ylim([-500, 500]);

% figure;
% subplot(2,1,1)
% plot(vel(1,:));
% subplot(2,1,2);
% plot(vel(2,:));

%% SAVE

% save("./sim.mat", 't', 'pos', 'vel', 'y', 'x_sensor');
save("./pf.mat", 't', 'pos', 'vel', 'y', 'x_sensor');

