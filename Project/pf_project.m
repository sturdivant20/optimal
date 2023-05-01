clc; clear; close all;
load("pf.mat");


n = 4;                      % number of states
m = size(y,1);              % number of measurements
s = size(x_sensor,1);       % number of sensors
L = size(t,2);              % number of time points
N = 1000;                   % number of particles

R = 3;                      % measurement noise [deg]
Q = diag([10, 10, 1, 1]);   % process noise [m, m, m/s, m/s]


% INITIALIZE
x = zeros(n,L);

% first update
p = [300.*randn(2,N); ...
     zeros(2,N)];
w = ones(1,N) ./ N;


% RUN SIMULATION
f = figure(Units='normalized', Position=[3.0, 0.4, 1.2, 0.4]);
for i = 1:L
    if i > 1
        dt = t(i) - t(i-1);
        [p, ~, w] = pf(p, w, x_sensor, y(:,i), Q, R, dt, i);
    end

    % plot
    plot(x_sensor(:,1), x_sensor(:,2), 'bo', LineWidth=3);
    hold on;
    plot(p(1,:), p(2,:), '.');
    for j = 1:2:size(pos,1)
        plot(pos(j,:), pos(j+1,:));
        plot(pos(j,1), pos(j+1,1), 'rx', LineWidth=3);
    end
    for j = 1:size(y,1)
        k = mod(j,size(x_sensor,1));
        if k == 0
            k = size(x_sensor,1);
        end
        plot([x_sensor(k,1), x_sensor(k,1) - 1000*cosd(y(j,i))], ...
             [x_sensor(k,2), x_sensor(k,2) - 1000*sind(y(j,i))], 'k');
    end
    hold off;
    xlim([-500, 500]);
    ylim([-500, 500]);
    legend('Sensor', 'Particles', 'Emitter', Location='northoutside', Orientation='horizontal');
    pause(0.1);
end

