function [xBar, Pxx] = ukf(x, P, y, x_sensor, Q, R, dt, qq)
% UKF unscented kalman filter
%
% Inputs:
%   x           4x1 current state (mean position and velocity)
%   P           4x4 current covariance
%   y           Lx1 angle of arrival measurements [deg]
%   x_sensor    Lx2 sensor state (position)
%   Q           LxL process noise
%   R           1x1 measurement noise [deg]
%   dt          1x1 delta time
%
% Outputs:
%   xBar        4x1 updated state mean
%   Pxx         4x4 error covariance
%
% Author:
%   Daniel Sturdivant   <dfs0012@auburn.edu>
%
% References:
%  [1]  E. A. Wan and R. Van Der Merwe, "The unscented Kalman filter for 
%       nonlinear estimation," Proceedings of the IEEE 2000 Adaptive 
%       Systems for Signal Processing, Communications, and Control 
%       Symposium (Cat. No.00EX373), Lake Louise, AB, Canada, 2000, pp. 
%       153-158, doi: 10.1109/ASSPCC.2000.882463. Online at: 
%       https://groups.seas.harvard.edu/courses/cs281/papers/unscented.pdf
%

m = size(y, 1);     % number of sensors
n = size(x, 1);     % number of states

y = deg2rad(y);     % measurement [rad]
R = deg2rad(R);     % measurement noise [rad]
x_old = x;          % previous state
xHat = zeros(n,2*n+1);
yHat = zeros(m,2*n+1);
xBar = zeros(n,1);
yBar = zeros(m,1);
Ux = zeros(m,2*n+1); 
Uy = zeros(m,2*n+1);
r = zeros(m,2*n+1);


% SCALING PARAMETERS
alpha = 1e-3;
beta = 2;
kappa = 0;
lambda = abs(alpha^2 * (n + kappa) - n);
Wm = zeros(1, 2*n+1);
Wc = zeros(1, 2*n+1);


% PROPAGATION
% Square root of covariance matrix
P_sqrt = chol((n + lambda) * P);

% state transition matrix (nonlinear w.r.t. time)
A = [1, 0, dt, 0; ...
     0, 1, 0, dt; ...
     0, 0, 1,  0; ...
     0, 0, 0,  1];

for i = 1:(2*n + 1)
    % propagate particles/ensembles with noise
    if i == 1
        xHat(:,i) = x_old;
        xHat(:,i) = A*xHat(:,i); % + Q*randn(n,1);
        Wm(i) = lambda / (n + lambda);
        Wc(i) = Wm(i) + (1 - alpha^2 + beta);
    elseif i > n+1
        xHat(:,i) = x_old - P_sqrt(:,i-n-1);
        xHat(:,i) = A*xHat(:,i); % + Q*randn(n,1);
        Wm(i) = 0.5 / (n + lambda);
        Wc(i) = Wm(i);
    else
        xHat(:,i) = x_old + P_sqrt(:,i-1);
        xHat(:,i) = A*xHat(:,i); % + Q*randn(n,1);
        Wm(i) = 0.5 / (n + lambda);
        Wc(i) = Wm(i);
    end

    % generate measurement from particle/ensemble
    dx = x_sensor - xHat(1:2,i)';
    Ux(:,i) = dx(:,1);
    Uy(:,i) = dx(:,2);
    r(:,i) = sum(dx.^2,2);
    yHat(:,i) =  wrapTo2Pi(atan2(Uy(:,i), Ux(:,i))); % + R*randn(m,1));

    % detect measurement ambiguity
    for j = 1:m
        if abs(y(j)-yHat(j,i)) > pi
            yHat(j,i) = mod(yHat(j,i),-2*pi);
        end
    end

    xBar = xBar + Wm(i) .* xHat(:,i);
    yBar = yBar + Wm(i) .* yHat(:,i);
end


Pxx = zeros(n,n);
Pyy = zeros(m,m);
Pxy = zeros(n,m);
for i = 1:(2*n + 1)
    Pxx = Pxx + Wc(i) .* (xHat(:,i) - xBar) * (xHat(:,i) - xBar)';   % state covariance
    Pyy = Pyy + Wc(i) .* (yHat(:,i) - yBar) * (yHat(:,i) - yBar)';   % measurement covariance
    Pxy = Pxy + Wc(i) .* (xHat(:,i) - xBar) * (yHat(:,i) - yBar)';   % state-measurement covariance
end
Pyy = Pyy + R;
Pxx = Pxx + Q;


% CORRECTION
L = Pxy * Pyy^-1;
Pxx = Pxx - L*Pyy*L';
xBar = xBar + L*(y - yBar);


end
