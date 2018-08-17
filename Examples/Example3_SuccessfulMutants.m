%% Example 3: Calculating the distribution of successful mutant arrival in two-type Belman-Harris Branching Process (BHBP)
% This example is continuation of the previous Example2. Here we discuss the joint probability that "successful mutant" 
% has not been born until time t and we do not have cells of type 1 (with subcritical reproduction rate, m_1 < 1) at 
% time t. This joint probability satisfies an integral equation that we will solve for different life length 
% distributions.
% P(T>t, Z^1(t)=0)=\int_0^t f_1(uq_0+(1-u)P(T>t-y, Z^1(t-y)=0)) dG_1(y).
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
T=100;      % defines the interval [0, T] in which we want to find numerical solution, i.e. the expected population count
h=0.1;      % defines the grid size of the numerical method

mean_life_length=10;    % we want to have average life length of 10
S=1-expcdf(0:h:T, mean_life_length)';   % defines the survivability function of the cells
H=[1-0.375, 0, 0.375]';      % defines the p.g.f of the offspring: P(0 offpsring) = 1-0.375, P(2 offspring) = 0.375, 
% or in other words the p.g.f. is 1-0.375 + 0*s^1 + 0.375*s^2. 

q_0=0.3;    % extinction probability
u=0.2;      % mutation probability

% defines the functions used in the integral equation
G=1-S;      % the life length distribution
z=zeros(size(S));
f_1=@(x)(arrayfun(@(x)(sum(H'.*x.^(0:size(H,1)-1))), x));   % the offspring p.g.f. for type 1 cells
f=@(y, Z)(f_1(u*q_0+(1-u)*Z));      % the equation as in article http://dx.doi.org/10.1016/j.csda.2016.12.013

Z_exp=numerical_solution(z, f, G, T, h);    % finds the solution

% calculates the solution for normally distributed life length
S=(1-normcdf(0:h:T, mean_life_length, 2.5)')./(1-normcdf(0, mean_life_length, 2.5));  % defines the survivability function of the cells
Z_normal=numerical_solution(z, f, 1-S, T, h);

% calculates the solution for uniformly distributed life length
S=1-unifcdf(0:h:T, 0, mean_life_length*2)';  % defines the survivability function of the cells
Z_unif=numerical_solution(z, f, 1-S, T, h);

% calculates the solution for beta distributed life length
S=[1-betacdf(0:h/(2*mean_life_length):1, 2, 2)'; zeros(T/h-length(0:h/(2*mean_life_length):1)+1, 1)];  % defines the survivability function of the cells
Z_beta=numerical_solution(z, f, 1-S, T, h);

% calculates the solution for chi^2 distributed life length
S=1-chi2cdf(0:h:T, mean_life_length)';  % defines the survivability function of the cells
Z_chi2=numerical_solution(z, f, 1-S, T, h);


%% make plots of the results
line_wd=2.5;
% plot the probability that successful mutant has not arrived yet
figure('visible','on', 'Units','pixels','OuterPosition',[0 0 1280 1024]);
set(gca,'FontSize',16)
hold on
plot(0:h:T, Z_exp, 'Color', [1, 0, 0, 0.5], 'LineWidth', line_wd)
plot(0:h:T, Z_normal, '-', 'Color', [0, 0, 0, 0.5], 'LineWidth', line_wd)
plot(0:h:T, Z_unif, '-', 'Color', [0, 1, 0, 0.5], 'LineWidth', line_wd)
plot(0:h:T, Z_beta, '-', 'Color', [0, 0, 1, 0.5], 'LineWidth', line_wd)
plot(0:h:T, Z_chi2, '--', 'Color', [1, 0, 0, 0.5], 'LineWidth', line_wd)
xlabel('Time'); ylabel('Probability P(T>t, Z^1(t)=0)');
legend('Exp(10) life length', 'N(10, 2.5) life length', 'U(0, 20) life length', 'Beta(2,2) life length, scaled in [0, 20]', ...
    'Chi^2(10) life length', 'Location', 'SE')
title('Successful mutant arrival in two-type Bellman-Harris BP')
print('./figures/Example3_fig1', '-dpng', '-r0')

