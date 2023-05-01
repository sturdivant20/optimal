function [xBar, Pxx] = ekfL(x, P, y, Q, R, dt, qq)
% EKFL loosely coupled extended kalman filter
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


% PROPAGATION
A = [1, 0, dt, 0; ...
     0, 1, 0, dt; ...
     0, 0, 1,  0; ...
     0, 0, 0,  1];

xBar = A*x;
Pxx = A*P*A' + Q;


% CORRECTION
C = [1, 0, 0, 0; ...
     0, 1, 0, 0];
L = Pxx*C'*(C*Pxx*C' + R)^-1;
Pxx = (eye(n) - L*C)*Pxx;
xBar = xBar + L*(y - C*xBar);

end
