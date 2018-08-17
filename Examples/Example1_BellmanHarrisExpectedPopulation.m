%% Example 1: Calculating the expected future population in Belman-Harris Branching Process (BHBP)
% The most used distribution for the life length in BHBP is exponential distribution, because it is the only case 
% in which we know the theoretical expectation for the population count. In all other cases we know that it satisfies 
% a certain renewal equation (Z(t) = S + \int_0^t [m*Z(t-y)] dG(y)) but we do not know its theoretical solution.
% This example considers the exponential life length and compares the numerical and theoretical solution. 
% Then it compares the results for exponential and truncated normal life length.
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
T=120;      % defines the interval [0, T] in which we want to find numerical solution, i.e. the expected population count
h=0.1;      % defines the grid size of the numerical method

mean_life_length=10;    % define the average life length to be 10
S=1-expcdf(0:h:T, mean_life_length)';   % defines the survivability function of the cells
H=[3/4,0,1/4]';     % defines the p.g.f of the offspring: P(0 offpsring) = 3/4, P(2 offspring) = 1/4, or in other words
% the p.g.f. is 3/4 + 0*s^1 + 1/4*s^2. 

% defines the functions used in the integral equation
G=1-S;      % the life length distribution
m=(0:(length(H)-1))*H;      % m=0.5 so the process is subcritical.
f=@(y, Z)(m*Z);     % the integrand in the equation
z=S;

Z_numerical_exp=numerical_solution(z, f, G, T, h);    % finds the solution, i.e. the expected population count in a Bellman-Harris branching process
Z_true=exp((m-1)*(0:h:T)/mean_life_length);     % the theoretical solution

% calculates the solution for normally distributed life length (truncated normal as the life length cannot be negative)
% here we do not have theoretical distribution, so numerical or simulation method is the only way to solve this
S=(1-normcdf(0:h:T, mean_life_length, 2.5)')./(1-normcdf(0, mean_life_length, 2.5));  % defines the survivability function of the cells
G=1-S;
z=S;
Z_numerical_normal=numerical_solution(z, f, G, T, h);

%% make plots of the results
% plot the numerical vs theoretical solution when the life length is exponentially distributed
line_wd=2.5;
figure('visible','on', 'Units','pixels','OuterPosition',[0 0 1280 1024]);
set(gca,'FontSize',16)
hold on
plot(0:h:T, Z_numerical_exp, 'Color', [1, 0, 0, 0.5], 'LineWidth', line_wd)
plot(0:h:T, Z_true, '--', 'Color', [0, 0, 0, 0.5], 'LineWidth', line_wd)
xlabel('Time'); ylabel('Expected population count');
legend('Numerical solution', 'Theoretical Solution')
title('Bellman-Harris BP with exponential life length Exp(10)')
print('./figures/Example1_fig1', '-dpng', '-r0')

% plot the expected population count fo exponentially distributed vs normally distributed life length
figure('visible','on', 'Units','pixels','OuterPosition',[0 0 1280 1024]);
set(gca,'FontSize',16)
hold on
plot(0:h:T, Z_numerical_exp, 'Color', [1, 0, 0, 0.5], 'LineWidth', line_wd)
plot(0:h:T, Z_numerical_normal, '--', 'Color', [0, 0, 0, 0.5], 'LineWidth', line_wd)
xlabel('Time'); ylabel('Expected population count');
legend('Exponential life length', 'Normal life length')
title('Bellman-Harris BP with mean life length of 10')
print('./figures/Example1_fig2', '-dpng', '-r0')


