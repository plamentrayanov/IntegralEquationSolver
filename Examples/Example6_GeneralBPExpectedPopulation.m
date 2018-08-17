%% Example 6: Solving Renewal Equations and estimating the renewal functions
% In this example is estimated the expected population size of the General Branching Process that has truncated normal
% distribution of the life length and normal shape of mu (mu is distribution function but not to a probability measure).
%
% It solves the equation Z(t) = S(t) + \int_0^t Z(t-y) d mu(y), where mu is the expectation of the point process, i.e. 
% a measure on [0, \infty], and S is the survivability function for females in the population, i.e. 1-S is the life
% length distribution.
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
addpath('../')  % adds the function numerical_solution to the path
T=250;          % defines the interval [0, T] in which we want to find numerical solution
h=0.01;         % defines the grid size of the numerical method

% Note: you can use the smoothed empirical distribution S, calculated from demographic data, instead of this truncated
% normal distribution. It is important however the last element to be 0
S=(1-normcdf(0:h:T,76, 10)')./(1-normcdf(0,76, 10));    % mean life length for women is 76
S(end)=0;   % survival probability for ages greater than omega is 0

m=0.7;      % the mean number of daughters a women has through her life
mu_women_pdf=normpdf(0:h:T, 28, 5)';
mu_women_pdf([1:12/h, 50/h:end])=0;     % women on ages less than 12 or greater than 50 do not give birth
mu=m*mu_women_pdf/(sum(mu_women_pdf)*h);    % scale mu so that it integrates to 0.7

z=S;
f=@(y, Z)(Z);   % defines the equation as renewal type equation

% the expected number of renewals satisfies this equation
Z_numerical=numerical_solution(z, f, cumsum(mu)*h, T, h);    % finds the solution, i.e. the expected population count in a Bellman-Harris branching process

% make plots of the results
line_wd=2.5;
figure('visible','on', 'Units','pixels','OuterPosition',[0 0 1280 1024]);
set(gca,'FontSize',16)
hold on
plot(0:h:T, Z_numerical, 'Color', [0, 0, 0.7, 0.5], 'LineWidth', line_wd)
xlabel('Time'); ylabel('Expected population count');
legend({'Numerical Solution to $Z(t) = S(t) + \int_0^t Z(t-y) d\mu(y)$'}, 'interpreter', 'latex')
title('Expected population count of Crump-Mode-Jagers BP')
print('./figures/Example6_fig1', '-dpng', '-r0')




