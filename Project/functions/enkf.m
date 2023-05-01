function [xHat, xBar, Pxx] = enkf(p, x, P, y, x_sensor, Q, R, dt, qq)
% ENKF ensemble kalman filter
%
% Inputs:
%   p           4xN particles state
%   x           4x1 current state (mean position and velocity)
%   P           4x4 current covariance
%   y           Lx1 angle of arrival measurements [deg]
%   x_sensor    Lx2 sensor state (position)
%   Q           LxL process noise
%   R           1x1 measurement noise [deg]
%   dt          1x1 delta time
%
% Outputs:
%   xHat        4xN updated particles state
%   xBar        4x1 updated state mean
%   Pxx         4x4 error covariance
%
% Author:
%   Daniel Sturdivant   <dfs0012@auburn.edu>
%
% References:
%  [1]  Matthias Katzfuss, Jonathan R. Stroud & Christopher K. Wikle (2016) 
%       Understanding the Ensemble Kalman Filter, The American 
%       Statistician, 70:4, 350-357, DOI: 10.1080/00031305.2016.1141709.
%       Online at: https://www.math.umd.edu/~slud/RITF17/enkf-tutorial.pdf
%

N = size(p, 2);     % number of particles
m = size(y, 1);     % number of sensors
n = size(x, 1);     % number of states

y = deg2rad(y);     % measurement [rad]
R = deg2rad(R);     % measurement noise [rad]
x_old = x;          % previous state
xHat = zeros(n,N);
yHat = zeros(m,N);
Ux = zeros(m,N); 
Uy = zeros(m,N);
r = zeros(m,N);


% PROPAGATION
% state transition matrix (nonlinear w.r.t. time)
A = [1, 0, dt, 0; ...
     0, 1, 0, dt; ...
     0, 0, 1,  0; ...
     0, 0, 0,  1];

for i = 1:N
    % propagate particles/ensembles with noise
    xHat(:,i) = A*p(:,i) + Q*randn(n,1);

    % generate measurement from particle/ensemble
    dx = x_sensor - xHat(1:2,i)';
    Ux(:,i) = dx(:,1);
    Uy(:,i) = dx(:,2);
    r(:,i) = sum(dx.^2,2);
    yHat(:,i) =  wrapTo2Pi(atan2(Uy(:,i), Ux(:,i))); % + R*randn);

    % detect measurement ambiguity
    for j = 1:m
        if abs(y(j)-yHat(j,i)) > pi
            yHat(j,i) = mod(yHat(j,i),-2*pi);
        end
    end
end

xBar = mean(xHat,2);
yBar = mean(yHat,2);
% Pxx = (xHat - xBar) * (xHat - xBar)';   % state covariance
Pyy = (yHat - yBar) * (yHat - yBar)';   % measurement covariance
Pyy = Pyy + R;
Pxy = (xHat - xBar) * (yHat - yBar)';   % state-measurement covariance


% CORRECTION
L = Pxy * Pyy^-1;
for i = 1:N
    xHat(:,i) = xHat(:,i) + L*(y - yHat(:,i));
end
xBar = mean(xHat,2);
Pxx = (xHat - xBar) * (xHat - xBar)';

% Pxx = Pxx - L*Pyy*L';
% x = x + L*(y - yHat);


end
