function [p, x, w] = pf(p, w, x_sensor, y, Q, R, dt, qq)
% PARTICLEFILTER nonlinear estimation using particle filter
%
% Inputs:
%   p           4xN particles state
%   x           4x1 current state (mean position and velocity)
%   w           1xN
%   x_sensor    Lx2 sensor state (position)
%   y           Lx1 angle of arrival measurements [deg]
%   Q           LxL process noise
%   R           1x1 measurement noise [deg]
%   dt          1x1 delta time
%
% Outputs:
%   x           4x1 updated state
%   p           4xM updated particles state
%   P           4x4 error covariance
%   
% Author:
%   Daniel Sturdivant   <dfs0012@auburn.edu>
%
% References:
%  [1]  Arulampalam et. al. (2002).  A tutorial on particle filters for 
%       online nonlinear/non-gaussian bayesian tracking. IEEE Transactions 
%       on Signal Processing. 50 (2). p 174--188. Online at: 
%       https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=978374
%
% Using Algorithm 3 from [1]
%


% INITIALIZE VARIABLES
L = size(y,1);      % number of sensors
N = size(p,2);      % number of particles
n = size(p,1);      % number of states
w_old = w;          % previous weights


% IMPORTANCE SAMPLING -> q(x[k] | x[k-1], y[k])
p_old = p;          % previous particles
p = zeros(n,N);     % new particles
w = zeros(1,N);     % new weights
yHat = zeros(L,N);

for i = 1:N
    % Propagation of state from PRIOR PDF (Eq. 62 [1]): 
    %   -> q(x[k] | x[k-1], y[k]) = p(x[k] | x[k-1]_i)
    %   -> probability of x given previous x
    A = [1, 0, dt, 0; ...
         0, 1, 0, dt; ...
         0, 0, 1,  0; ...
         0, 0, 0,  1];
    p(:,i) = A*p_old(:,i) + Q*randn(4,1);

    % Weights from from prior PDF (Eq. 63 [1])
    %   -> w[k] = w[k-1] * p(y[k] | x[k])
    for j = 1:L
        k = mod(j,size(x_sensor,1));
        if k == 0
            k = size(x_sensor,1);
        end
        dx = x_sensor(k,:) - p(1:2,i)';
        yHat(j,i) = wrapTo360(atan2d(dx(2), dx(1))); % azimuth observation to j-th sensor
        if abs(y(j) - yHat(j,i)) > 180
            yHat(j,i) = mod(yHat(j,i),-360);
        end
        w(i) = w_old(i) * normpdf(y(j) - yHat(j,i), 0, R);  % normal distribution
    end
end

% Normalize weights and calculate effective sample size
w = w ./ sum(w);
N_eff = 1 / sum(w.^2);  % (Eq. 51 [1])
threshold = 0.5 * N;
if N_eff < threshold
    [p, w, ~] = resample(p, w, N);
end

% Ouptut state, particles, and weight
x = sum(w.*p,2);
P = zeros(n,n);


% RESAMPLING -> Algorithm 2 from [1]
function [x, w, idx] = resample(x, w, N)

% Generate CDF
C = zeros(size(w));
for ii = 2:N
    C(ii) = C(ii-1) + w(ii);
end

ii = 1;
idx = zeros(1,N);   % parent indexes
u1 = rand/N;        % draw a starting point
for jj = 1:N
    uj = u1 + 1/N * (jj-1);  % move along CDF
    while uj > C(ii) && (ii < N)
        ii = ii+1;
    end
    x(:,jj) = x(:,ii);  % assign sample
    w(jj) = w(ii);      % assign weight
    idx(jj) = ii;       % assign parent
end

w = w ./ sum(w);

end

end