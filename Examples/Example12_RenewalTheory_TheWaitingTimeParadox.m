%% Example 12: The waiting time paradox. Solving the Renewal Equations for age, residual lifetime and total lifetime
% This example replicates the results in Example 11, but for truncated normal distribution of inter-arrival times.
%
%
%       Copyright (C) 2018  Plamen Ivaylov Trayanov
%
%       This program is free software: you can redistribute it and/or modify
%       it under the terms of the GNU General Public License as published by
%       the Free Software Foundation, either version 3 of the License, or
%       (at your option) any later version.
% 
%       This program is distributed in the hope that it will be useful,
%       but WITHOUT ANY WARRANTY; without even the implied warranty of
%       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%       GNU General Public License for more details.
% 
%       You should have received a copy of the GNU General Public License
%       along with this program. If not, see <http://www.gnu.org/licenses/>.
%
%
addpath('../')      % adds the function numerical_solution to the path
T=50;          % defines the interval [0, T] in which we want to find numerical solution, i.e. the expected population count
h=0.05;         % defines the grid size of the numerical method

%% define the renewal process with normal lifetime
mean_life_length=10;    % we want to have an average renewal time of 10
sigma_life_length=2.5;
F=normcdf(0:h:T, mean_life_length, sigma_life_length)'./(1-normcdf(0, mean_life_length, sigma_life_length));   % defines the distributions of the time between renewals to be normal
f=@(y, Z)(Z);   % defines the equation as renewal-type equation

T_A=floor(T/2);         % the maximum possible horizon for evaluating F_B (see the formula)
F_A=zeros(T_A/h+1);
F_B=zeros(size(F_A));
F_C=zeros(size(F_A));
for i=1:size(F_A,2)
    % calculate P(A(t)<=x(i))
    z=(1-F(1:T_A/h+1)) .* [ones(i,1); zeros(T_A/h-i+1, 1)];
    F_A(:,i)=numerical_solution(z, f, F(1:T_A/h+1), T_A, h);
    
    % calculate P(B(t)>x)
    z=(1-F(i:(i+T_A/h)));
    F_B(:,i)=1-numerical_solution(z, f, F(1:T_A/h+1), T_A, h);
    
    % calculate P(C(t)<=x)
    z=(F(i) - F(1:T_A/h+1)) .* [ones(i,1); zeros(T_A/h-i+1, 1)];
    F_C(:,i)=numerical_solution(z, f, F(1:T_A/h+1), T_A, h);
end
dF_A=diff(F_A, 1, 2)/h;
dF_B=diff(F_B, 1, 2)/h;
dF_C=diff(F_C, 1, 2)/h;

%% calculate the theoretical limit distributions
mu=sum(1-F)*h;             % mu = E(T) = \int_0^\infty 1-F(u) du
F_0=cumsum(1-F)*h/mu;      % F_0(x) = 1/mu * \int_0^x 1-F(u) du
F_A_lim=F_0;
F_B_lim=F_0;
F_C_lim=cumsum((0:h:T-h)' .* diff(F))/mu;      % 1/mu * \int_0^x u dF(u)

dF_A_lim=diff(F_A_lim)/h;
dF_B_lim=diff(F_B_lim)/h;
dF_C_lim=diff(F_C_lim)/h;

%% make plots of A(t), B(t) and C(t) for exponential inter-arrival time distribution
line_wd=2.5;
figure('visible','on', 'Units','pixels','OuterPosition',[0 0 1600 800]);
set(gca,'FontSize',16)
colormap(jet)
hold on
subplot(1,2,1), mesh(h:10*h:T_A, 0:10*h:T_A, dF_A(1:10:end, 1:10:end), 'FaceAlpha',1, 'FaceColor', 'interp');
material('shiny'), lighting('gouraud'), light, view(3)
xlabel('x'); ylabel('t'); zlabel('$P(A(t) \in [x, x+dx))$','interpreter','latex')
title(strcat('Renewal Process with N(', num2str(mean_life_length), ', ', num2str(sigma_life_length), ')', ...
    ' arrival time. P.D.F of A(t) - Age'),'interpreter','latex');
subplot(1,2,2), mesh(0:10*h:T_A, 0:10*h:T_A, F_A(1:10:end, 1:10:end), 'FaceAlpha',1, 'FaceColor', 'interp');
material('shiny'), lighting('gouraud'), light, view(3)
xlabel('x'); ylabel('t'); zlabel('$P(A(t)\leq x)$','interpreter','latex')
title({'C.D.F of $A(t)$: $P(A(t)\leq x) = (1-F(t))1_{[0, x]}(t) + \int \limits_0^t P(A(t-y)\leq x) \, \mathrm{d} F(y)$'},'interpreter','latex');
print('./figures/Example12_fig1', '-dpng', '-r0')

figure('visible','on', 'Units','pixels','OuterPosition',[0 0 1600 800]);
set(gca,'FontSize',16)
colormap(jet)
hold on
subplot(1,2,1), mesh(h:10*h:T_A, 0:10*h:T_A, dF_B(1:10:end, 1:10:end), 'FaceAlpha',1, 'FaceColor', 'interp');
material('shiny'), lighting('gouraud'), light, view(3)
xlabel('x'); ylabel('t'); zlabel('$P(B(t) \in [x, x+dx))$','interpreter','latex')
title(strcat('Renewal Process with N(', num2str(mean_life_length), ', ', num2str(sigma_life_length), ')', ...
    ' arrival time. P.D.F of B(t) - residual lifetime'),'interpreter','latex');
subplot(1,2,2), mesh(0:10*h:T_A, 0:10*h:T_A, F_B(1:10:end, 1:10:end), 'FaceAlpha',1, 'FaceColor', 'interp');
material('shiny'), lighting('gouraud'), light, view(3)
xlabel('x'); ylabel('t'); zlabel('$P(B(t) > x)$','interpreter','latex')
title({'C.D.F of $B(t)$: $P(B(t) > x) = (1-F(t+x)) + \int \limits_0^t P(B(t-y) > x) \, \mathrm{d} F(y)$'},'interpreter','latex');
print('./figures/Example12_fig2', '-dpng', '-r0')

figure('visible','on', 'Units','pixels','OuterPosition',[0 0 1600 800]);
set(gca,'FontSize',16)
colormap(jet)
hold on
subplot(1,2,1), mesh(h:5*h:T_A, 0:5*h:T_A, dF_C(1:5:end, 1:5:end), 'FaceAlpha',1, 'FaceColor', 'interp');
material('shiny'), lighting('gouraud'), light, view(3)
xlabel('x'); ylabel('t'); zlabel('$P(C(t) \in [x, x+dx))$','interpreter','latex')
title(strcat('Renewal Process with N(', num2str(mean_life_length), ', ', num2str(sigma_life_length), ')', ...
    ' arrival time. P.D.F. of C(t) - total lifetime'),'interpreter','latex');
subplot(1,2,2), mesh(0:5*h:T_A, 0:5*h:T_A, F_C(1:5:end, 1:5:end), 'FaceAlpha',1, 'FaceColor', 'interp');
material('shiny'), lighting('gouraud'), light, view(3)
xlabel('x'); ylabel('t'); zlabel('$P(C(t) \leq x)$','interpreter','latex')
title({'C.D.F of $C(t)$: $P(C(t) \leq x) = (F(x)-F(t))1_{[0, x]}(t) + \int \limits_0^t P(C(t-y) \leq x) \, \mathrm{d} F(y)$'},'interpreter','latex');
print('./figures/Example12_fig3', '-dpng', '-r0')

%% plot the waiting time paradox and the limiting distribution
figure('visible','on', 'Units','pixels','OuterPosition',[0 0 1280 1024]);
set(gca,'FontSize',16)
hold on
plot(h:h:T_A, dF_C(end,:), 'Color', [0.7, 0, 0, 0.7], 'LineWidth', line_wd);
plot(h:h:T_A, dF_C_lim(1:T_A/h), '--', 'Color', [0, 0, 0, 0.7], 'LineWidth', 3);
plot(h:h:T_A, diff(F(1:T_A/h+1))/h, '--', 'Color', [0, 0, 0.7, 0.7], 'LineWidth', line_wd);
xlabel('x'); ylabel('Probability')
legend({'$P(C(t) \in [x, x+\mathrm{d}x))$, for t=25', '$\lim \limits_{t \rightarrow \infty} P\left(C(t) \in [x, x+\mathrm{d}x)\right)=\frac{\mathrm{d}}{\mathrm{d}x}\left(\frac{1}{\mu} \int \limits_0^x u \, \mathrm{d}F(u)\right)$', 'Inter-arrival time distribution, N(10, 2.5)'},'interpreter','latex', 'Location', 'NE')
text(15, 0.12, ['The distribution of C(t) is different than \newline', ...
    'the inter-arrival time distribution.'], 'FontSize', 14)
print('./figures/Example12_fig4', '-dpng', '-r0')

mean_B=sum(1-F_B(end,:))*h;
mean_F=sum(1-F)*h;
figure('visible','on', 'Units','pixels','OuterPosition',[0 0 1280 1024]);
set(gca,'FontSize',16)
hold on
plot(h:h:T_A, dF_B(end,:), 'Color', [0.7, 0, 0, 0.7], 'LineWidth', line_wd);
plot(h:h:T_A, dF_B_lim(1:T_A/h), '--', 'Color', [0, 0, 0, 0.7], 'LineWidth', 3);
plot(h:h:T_A, diff(F(1:T_A/h+1))/h, '--', 'Color', [0, 0, 0.7, 0.7], 'LineWidth', line_wd);
xlabel('x'); ylabel('Probability')
legend({'$P(B(t) \in [x, x+\mathrm{d}x))$, for t=25', '$\lim \limits_{t \rightarrow \infty} P(B(t) \in [x, x+\mathrm{d}x))=\frac{\mathrm{d}}{\mathrm{d}x}\left(\frac{1}{\mu} \int \limits_0^x 1-F(u) \, \mathrm{d}u\right)$', 'Distribution of inter-arrival time, N(10, 2.5)'},'interpreter','latex', 'Location', 'NE')
text(10, 0.06, ['The paradox: The average waiting time to the next renewal is ', num2str(mean_B), ',', newline,...
    'greater than half of the average inter-renewal time (', num2str(mean_F/2),').', newline, ...
    'For large t, the theoretical limit is $(10^2+2.5^2)/(2 \cdot 10) = 5.31$.'],...
    'FontSize', 14, 'interpreter', 'latex')
title('The waiting time paradox for N(10, 2.5) inter-arrival time.')
print('./figures/Example12_fig5', '-dpng', '-r0')




