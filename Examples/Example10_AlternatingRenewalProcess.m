%% Example 10: Solving the renewal equation for Alternating Renewal Process and estimating the renewal function
% The example follows the lecture by Professor Whitt published in http://www.columbia.edu/~ww2040/3106F13/lect1119.pdf
% IEOR 3106: Fall 2013, Professor Whitt, Alternating Renewal Processes and The Renewal Equation 
% 
% We have a computer working on 3 parts that can brake and take time to repair. The computer stops working if any part
% fails. During repair, the other parts do not wear.
%
% U1, U2, U3 - time before failing for each part
% D1, D2, D3 - time needed to repair each part
%
% Z(t) - state at time t: up or down
% In the lecture example U1 ~ Exp(10), U2 ~ Exp(20) and U3 ~ Exp(30). 
% U = min(U1, U2, U3) - total uptime is Exp(1/(1/10+1/20+1/30)), i.e. Exp(60/11)
%
% D1 - Exp(1) lifetime
% D2 - Unif(1,2) lifetime
% D3 - Gamma(k, theta) lifetime, where we know E=k*theta=3, STD=sqrt(k)*theta=10, so k=0.09, theta=33.3333
%
% In the lecture, there are calculated probabilities for each of the 3 parts to fail first:
% P(N=1)=6/11, P(N=2)=3/11 and P(N=3)=6/11.
%
% D - F_D(t) = P(D<=t) = P(D<=t | N=1) * P(N=1) + P(D<=t | N=1) * P(N=2) + P(D<=t | N=1) * P(N=3)=
%            = 6/11 * Exp(1) + 3/11 * Unif(1,2) + 2/11 * Gamma(k, theta)
%
% F_U(t) denotes the CDF of U, i.e. Exp(60/11)
% F(t) denotes the CDF of U + D, which is a convolution of F_U * F_D
%
% The the probability for the system to be up at time t satisfies the renewal equation:
% P(Z(t)=1) = 1-F_U(t) + \int_0^t P(Z(t-s)=1) dF(s)
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
T=100;              % defines the interval [0, T] in which we want to find numerical solution, i.e. the expected population count
h=0.01;             % defines the grid size of the numerical method

F_U=expcdf(0:h:T, 1/(1/10 + 1/20 + 1/30))';   % distribution of the uptime U

% P(D <=t) = P(D<=t | N=1) * P(N=1) + P(D<=t | N=1) * P(N=2) + P(D<=t | N=1) * P(N=3)=
% = 6/11 * Exp(1) + 3/11 * Unif(1,2) + 2/11 * Gamma(k, theta)
F_D = 6/11 * expcdf(0:h:T, 1)' + 3/11 * unifcdf(0:h:T, 1, 3)' + 2/11 * gamcdf(0:h:T, 0.09, 33.3333)';

F=convolution(F_U, F_D);   % distribution of U + D

f=@(y, Z)(Z);   % defines the equation as renewal-type equation
P_numerical=numerical_solution(1-F_U, f, F, T, h);

% the expected number of renewals, including the initial moment
m_t=numerical_solution(ones(T/h+1,1), f, F, T, h);  % estimates the renewal function as a solution to renewal equation

%% plot the distributions of Uptime, Downtime and Inter-arrival time (Uptime+Downtime)
line_wd=2.5;
figure('visible','on', 'Units','pixels','OuterPosition',[0 0 1280 1024]);
set(gca,'FontSize',16)
hold on
plot(h:h:T, diff(F_U)/h, 'Color', [0, 0, 0.7 0.6], 'LineWidth', line_wd)
plot(h:h:T, diff(F_D)/h, 'Color', [0.7, 0, 0 0.6], 'LineWidth', line_wd)
plot(h:h:T, diff(F)/h, 'Color', [0, 0, 0 0.6], 'LineWidth', line_wd)
xlim([0, 10])
ylim([0, 1])
xlabel('Time'); ylabel('Probability');
legend('P.D.F. of U', 'P.D.F. of D', 'P.D.F. of U + D')
title('Alternating Renewal Process')
print('./figures/Example10_fig1', '-dpng', '-r0')

%% plot the probability the system is up at time t
figure('visible','on', 'Units','pixels','OuterPosition',[0 0 1280 1024]);
set(gca,'FontSize',16)
hold on
plot(0:h:T, P_numerical, 'Color', [0, 0, 0.7 0.7], 'LineWidth', line_wd)
plot([0 T], [60/78, 60/78], 'r', 'LineWidth', line_wd)    % E(U)/E(U+D) = 60/78, calculated in the lecture example
xlabel('Time'); ylabel('Probability');
legend('Probability the system is up', 'Theoretical limit')
title('Alternating Renewal Process')
print('./figures/Example10_fig2', '-dpng', '-r0')

%% plot the probability the system is up at time t
figure('visible','on', 'Units','pixels','OuterPosition',[0 0 1280 1024]);
set(gcf, 'defaultAxesColorOrder',[0, 0, 0.7; 0.7, 0, 0])
set(gca, 'FontSize', 16)
hold on
yyaxis left
plot(0:h:T, m_t, 'Color', [0, 0, 0.7 0.7], 'LineWidth', line_wd)
ylabel('m(t)');
yyaxis right
plot(h:h:T, diff(m_t)/h, 'Color', [0.7, 0, 0, 0.7], 'LineWidth', line_wd)
ylabel('m''(t)');
xlabel('Time'); 
legend('Expected number of renewals', 'Derivative')
title('Alternating Renewal Process')
print('./figures/Example10_fig3', '-dpng', '-r0')

% m_t resembles a linear function but it is not linear, as seen from the derivative

%% private functions to this script
function res=convolution(F1, F2)
% calculates the convolution of two functions presented as vectors
res=zeros(size(F1));
if iscolumn(F1)
    dF2=[diff(F2); F2(end)-F2(end-1)];
else
    dF2=[diff(F2), F2(end)-F2(end-1)];
end
for i=1:length(res)
    res(i) = sum(F1(i:-1:1) .* dF2(1:i));
end
end

