%% Daniel Sturdivant | Optimal HW1
clc; clear; close all;

f = figure('units','normalized','position',[0.1 0.1 0.3 0.7]);
% f = figure('units','normalized','position',[0.1 0.1 0.4 0.5]);
tabs = uitabgroup(f);
tab(1) = uitab(Title='1)');
tab(2) = uitab(Title="2)");
tab(3) = uitab(Title="3)");
tab(4) = uitab(Title="4)");
tab(5) = uitab(Title="5)");
tab(6) = uitab(Title="6)");
tab(7) = uitab(Title="7)");


%% Problem 1
fprintf("\n<strong>PROBLEM 1</strong>\n");

% probability mass functions
prob1.pmf_a = [1/6, 1/6, 1/6, 1/6, 1/6, 1/6];
prob1.pmf_b = [1/6, 1/6, 1/6, 1/6, 1/6, 1/6];
prob1.pmf_c = [1/3, 0, 1/2, 0, 1/6, 0];
prob1.pmf_d_1 = [1/6, 1/6, 1/6, 1/6, 1/6, 1/6];
prob1.pmf_d_2 = [1/3, 0, 1/2, 0, 1/6, 0];

% initialize
prob1.a = prob1.pmf_a;
prob1.b = prob1.pmf_b;
prob1.c = prob1.pmf_c;
prob1.d = prob1.pmf_d_1;

% convolve pdf for 5 more throws
for i = 1:5
    prob1.a = conv(prob1.a, prob1.pmf_a);
    prob1.b = conv(prob1.b, prob1.pmf_b);
    prob1.c = conv(prob1.c, prob1.pmf_c);
    if (i < 3)
        prob1.d = conv(prob1.d, prob1.pmf_d_1);
    else
        prob1.d = conv(prob1.d, prob1.pmf_d_2);
    end
end

% check sum == 1
prob1.sum_a = sum(prob1.a);
prob1.sum_b = sum(prob1.b);
prob1.sum_c = sum(prob1.c);
prob1.sum_d = sum(prob1.d);

prob1.ax = axes(tab(1));
hold(prob1.ax, "on");
stem(prob1.ax, 6:36, prob1.a*100, LineWidth=3, DisplayName="Part A");
stem(prob1.ax, 24:54, prob1.b*100, LineWidth=3, DisplayName="Part B");
stem(prob1.ax, 6:36, prob1.c*100, LineWidth=3, DisplayName="Part C");
stem(prob1.ax, 6:36, prob1.d*100, LineWidth=3, DisplayName="Part D"); 
grid(prob1.ax, "on");
legend(prob1.ax);
title(prob1.ax, "\textbf{PDFs of 6 Dice Rolls}", Interpreter="latex");
xlabel(prob1.ax, "Sum of Dice Rolls");
ylabel(prob1.ax, "Probability (%)");
prob1.ax.YAxis.TickLabelFormat = "%g%%";
prob1.ax.FontSize = 16;

% exportgraphics(tab(1), './media/p1.png');

%% Problem 2
fprintf("\n<strong>PROBLEM 2</strong>\n");

% D1 = D2, f_D1 = f_D2
prob2.D1   = [ 1 ,  2 ,  3 ,  4 ,  5 ,  6 ];
prob2.D1_f = [1/6, 1/6, 1/6, 1/6, 1/6, 1/6];
prob2.D1D2_f = prob2.D1_f' * prob2.D1_f;

% part a
prob2.D1_mu = 1/length(prob2.D1) * sum(prob2.D1);
prob2.D1_2 = 1/length(prob2.D1) * sum(prob2.D1.^2);
prob2.D1_var = prob2.D1_2 - prob2.D1_mu^2;

% part b
prob2.D1D2_cov = [prob2.D1_var,0;0,prob2.D1_var];

% part c
prob2.v1 = [1,2,3,4,5,6];
prob2.v2 = [2,3,4,5,6,7,8,9,10,11,12];
prob2.v1_f = [1/6, 1/6, 1/6, 1/6, 1/6, 1/6];
prob2.v2_f = [1/36, 2/36, 3/36, 4/36, 5/36, 6/36, 5/36, 4/36, 3/36, 2/36, 1/36];
prob2.v1v2_f = prob2.v1_f' * prob2.v2_f;

% part d
prob2.v1_mu = 1/length(prob2.v1) * sum(prob2.v1);
prob2.v1_2 = 1/length(prob2.v1) * sum(prob2.v1.^2);
prob2.v1_var = prob2.v1_2 - prob2.v1_mu^2;

% part e
prob2.v2_mu = sum(prob2.v2_f.*prob2.v2);
prob2.v2_2 = sum(prob2.v2_f.*(prob2.v2.^2));
prob2.v2_var = prob2.v2_2 - prob2.v2_mu^2;

% part f
prob2.v1v2_cov = [prob2.v1_var, prob2.v1_var; prob2.v1_var, prob2.v1_var+prob2.v2_var];


%% Problem 3
fprintf("\n<strong>PROBLEM 3</strong>\n");


%% Problem 4
fprintf("\n<strong>PROBLEM 4</strong>\n");

prob4.dice = [-2.5, -1.5, -0.5, 0.5, 1.5, 2.5];
prob4.sums = -5:5;

% part a
prob4.v0_f = conv([1/6, 1/6, 1/6, 1/6, 1/6, 1/6], [1/6, 1/6, 1/6, 1/6, 1/6, 1/6]);

% part b
prob4.v0_mu = sum(prob4.v0_f .* prob4.sums);
prob4.v0_var = sum(prob4.v0_f .* prob4.sums.^2);

% part c -> vn[k+1] = (1-r)*vn[k] + r*v0[k]
i = 1;
prob4.VN_var = zeros(100,1);
prob4.VN_mean = zeros(100,1);
for r = linspace(0,1.9,100)
    prob4.V0 = zeros(10000,1);
    prob4.VN = zeros(10000,1);
    for k = 1:10000
        prob4.V0(k) = sum([prob4.dice(randi(6)), prob4.dice(randi(6))]);
        prob4.VN(k+1) = (1-r)*prob4.VN(k) + r*prob4.V0(k);
    end
    prob4.VN_var(i) = var(prob4.VN);
    prob4.VN_mu(i) = mean(prob4.VN);
    i = i + 1;
end

R = zeros(100,1);
R(1) = 1;
r = 0.1;
for i = 1:100
    R(i+1) = (1-r)^2 * R(i);
end

prob4.ax1 = axes(Parent=tab(3));
plot(prob4.ax1, 0:100, R, LineWidth=2.5);
grid(prob4.ax1, "on");
title(prob4.ax1, "\textbf{Autocorrelation}", Interpreter="latex");
xlabel(prob4.ax1, "R(k)");
ylabel(prob4.ax1, "Number of Dice Rolls");
prob4.ax1.FontSize = 16;

exportgraphics(tab(3), './media/p4-1.png');

prob4.ax = axes(Parent=tab(4));
hold(prob4.ax, "on");
plot(prob4.ax, linspace(0,1.9,100), prob4.VN_var, '-', LineWidth=2.5, DisplayName='Variance');
plot(prob4.ax, linspace(0,1.9,100), prob4.VN_mu, '-', LineWidth=2.5, DisplayName='Mean');
legend(prob4.ax, Location="northwest");
grid(prob4.ax, "on");
title(prob4.ax, "\textbf{$V_N$ Statistics Given $r$}", Interpreter="latex");
xlabel(prob4.ax, "r");
ylabel(prob4.ax, "Mean/Variance");
prob4.ax.FontSize = 16;

% exportgraphics(tab(4), './media/p4.png');


%% Problem 5
fprintf("\n<strong>PROBLEM 5</strong>\n");

prob5.f1 = @(x) x.*x./2;
prob5.mu = integral(prob5.f1,0,2);

prob5.f2 = @(x) x.^2.*x./2;
prob5.var = integral(prob5.f2,0,2) - prob5.mu^2;


%% Problem 6
fprintf("\n<strong>PROBLEM 6</strong>\n");

% zero mean
prob6.P = [2, 1; 1, 4];

% part a
[prob6.v, prob6.s] = eig(prob6.P);

% part b
prob6.A = prob6.s^(-1/2) * prob6.v;

[~, part6.idx] = max(diag(prob6.s));
prob6.alpha = atand(prob6.v(1,part6.idx) / prob6.v(2,part6.idx)); % rotation angle (ccw)
prob6.RotMat = [cosd(prob6.alpha), sind(prob6.alpha); -sind(prob6.alpha), cosd(prob6.alpha)];
prob6.principle_axes = diag(prob6.s^(1/2)); % principle axes

% part c
prob6.theta = 0:360;
prob6.a = [cosd(prob6.theta); sind(prob6.theta)];
prob6.b1 = prob6.A \ (0.25 .* prob6.a);
prob6.b2 = prob6.A \ (1.00 .* prob6.a);
prob6.b3 = prob6.A \ (1.50 .* prob6.a);

prob6.ax = axes(Parent=tab(6));
hold(prob6.ax, "on");
fill(prob6.ax, prob6.b3(1,:), prob6.b3(2,:), [0.64 0.08 0.18], FaceAlpha=0.3, EdgeColor=[0.64 0.08 0.18], LineWidth=2.5, DisplayName="k=1.50");
fill(prob6.ax, prob6.b2(1,:), prob6.b2(2,:), [0.00 0.45 0.74], FaceAlpha=0.3, EdgeColor=[0.00 0.45 0.74], LineWidth=2.5, DisplayName="k=1.00");
fill(prob6.ax, prob6.b1(1,:), prob6.b1(2,:), [0.47 0.67 0.19], FaceAlpha=0.3, EdgeColor=[0.47 0.67 0.19], LineWidth=2.5, DisplayName="k=0.25");
quiver(prob6.ax, 0, 0, prob6.v(1,1), prob6.v(2,1), 'k', LineWidth=2.5, HandleVisibility="off");
quiver(prob6.ax, 0, 0, prob6.v(1,2), prob6.v(2,2), 'k', LineWidth=2.5, HandleVisibility="off");
grid(prob6.ax, "on");
legend(prob6.ax, Location="northwest");
title(prob6.ax, "\textbf{Error Ellipses}", Interpreter="latex");
xlabel(prob6.ax, "x");
ylabel(prob6.ax, "y");
prob6.ax.FontSize = 16;

% exportgraphics(tab(6), './media/p6.png');

% part d (linear interpolation from error ellipse document)
prob6.p1 = 0 + (0.25-0)*(39.34-0)/(1-0);
prob6.p2 = 39.34;
prob6.p3 = 50 + (1.5-1.177)*(90-50)/(2.146-1.177);


%% Problem 7
fprintf("\n<strong>PROBLEM 7</strong>\n");

prob7.x = linspace(-8,8,1000);  % x values
prob7.y = linspace(-8,8,1000);    % x values
prob7.sigma = 2;                % x standard deviation
prob6.mu = 0;                   % x mean

prob7.pdf_x = (1/sqrt(2*pi*prob7.sigma^2)) .* exp(-prob7.x.^2 ./ (2*prob7.sigma^2));
prob7.pdf_y = abs(real(1./(2 .* sqrt(prob7.y./2)))) .* (1/sqrt(2*pi*prob7.sigma^2)) ...
                .* exp(-prob7.y ./ (4*prob7.sigma^2));

prob7.ax = axes(Parent=tab(7));
hold(prob7.ax, "on");
plot(prob7.x, prob7.pdf_x*100, LineWidth=2.5, DisplayName="f_X");
plot(prob7.y, prob7.pdf_y*100, LineWidth=2.5, DisplayName="f_Y");
grid(prob7.ax, "on");
legend(prob7.ax);
title(prob7.ax, "\textbf{PDF's for X and Y}", Interpreter="latex");
ylabel(prob7.ax, "Probability (%)");
xlabel(prob7.ax, "Value");
xlim(prob7.ax, [-8 8]);
prob7.ax.YAxis.TickLabelFormat = "%g%%";
prob7.ax.FontSize = 16;

% exportgraphics(tab(7), './media/p7.png');



