function [x, P] = ls(x_sensor, y, x, R)

L = length(y);
y = deg2rad(y);
i = 0;
error = 1 ;
lambda = 1e-5;

while error > 1e-6 && i < 100
    
    u = x_sensor - x(1:2)';
    r = sum(u.^2,2);
    yHat = wrapTo2Pi(atan2(u(:,2), u(:,1)));

    H = [u(:,2), -u(:,1)] ./ r;
    for j = 1:L
        if abs(y(j)-yHat(j)) > pi
            yHat(j) = mod(yHat(j),-2*pi);
        end
    end
    dy = y - yHat;

    % Levenburg Update
    dx = (H'*H + lambda*eye(2))^-1*H'*dy;
    tmp = sqrt(sum(dx.^2,"all"));
    if tmp < error
        lambda = lambda / 10;
        x(1:2) = x(1:2) + dx;
        error = tmp;
    else
        lambda = lambda * 10;
    end

    i = i+1;
end

P = deg2rad(R) .* (H'*H)^-1;

end