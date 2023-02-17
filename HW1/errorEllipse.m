function [b50, b90, b99] = errorEllipse(P, mu)
% ERRORELLIPSE creates 2D error ellipse given 2x2 covariance matrix
%
% Inputs:
%   P       2x2 Covariance matrix between X and Y
%   mu      2x1 Mean of X and Y(optional)
%
% Outputs:
%   b50     2x361 X and Y of error ellipse (50% confidence region)
%   b90     2x361 X and Y of error ellipse (50% confidence region)
%   b99     2x361 X and Y of error ellipse (50% confidence region)
%
% Author:
%   Daniel Sturdivant   <dfs0012@auburn.edu>
%

if nargin < 2
    mu = [0;0];
end

% 1) find the eigenvalues and eigenvectors of P
[V,D] = eig(P);

% 2) Compute the transform matrix
Ainv = inv(D^(-1/2) * V);

% 3) Generate Circle
theta = 0:360;
a = [cosd(theta); sind(theta)];

% 4) Transform into ellipse of deisred radius
b50 = mu + (Ainv * (1.177 * a));
b90 = mu + (Ainv * (2.146 * a));
b99 = mu + (Ainv * (3.035 * a));

end