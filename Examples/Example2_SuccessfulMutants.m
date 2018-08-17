%% Example 2: Calculating the distribution of successful mutant arrival in two-type Belman-Harris Branching Process (BHBP)
% This example considers the model in "Branching processes in continuous time as models of mutations: Computational 
% approaches and algorithms by Slavtchova-Bojkova, M., Trayanov, P., Dimitrov, S. (2017). Computational Statistics and 
% Data Analysis, http://dx.doi.org/10.1016/j.csda.2016.12.013.
%
% Short description of the article model:
% type 0 particles - supercritical, do not mutate to other types once they occur
% type 1 particles - subcritical, with probability of mutation to type 0
% The population starts from 1 particle of type 1.
% "Successful mutant" - a particle of type 0 that starts a type 0 branching, that never extincts
% "Unsuccessful mutant" - a particle of type 0 that starts a type 0 branching, that extincts
% The object of interest in the article is to calculate the distribution T of the first arrival of a "successful mutant".
%
% The distribution of the arrival time satisfies the integral equation:
% P(T>t) = 1-G(t) + \int_0^t f_1(u*q_0+(1-u)*P(T>t-y)) dG(y)
% The example shows the results for different distributions of life lengths.
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
H=[1-0.375, 0, 0.375]';      % defines the p.g.f of the offspring: P(0 offpsring) = 0.6250, P(2 offspring) = 0.375, or in other words
% the p.g.f. is 0.6250 + 0*s^1 + 0.375*s^2. 

q_0=0.3;    % extinction probability for type 0
u=0.2;      % mutation probability from type 1 to type 0

% defines the functions used in the integral equation
G=1-S;      % the life length distribution
z=S;
f_1=@(x)(arrayfun(@(x)(sum(H'.*x.^(0:size(H,1)-1))), x));   % the offspring p.g.f. for type 1 cells
f=@(y, Z)(f_1(u*q_0+(1-u)*Z));      % defines the integrand in the integral equation

Z_exp=numerical_solution(z, f, G, T, h);    % finds the solution, i.e. P(T>t), for exponential distribution G(t)

% calculates the solution for normally distributed life length
S=(1-normcdf(0:h:T, mean_life_length, 2.5)')./(1-normcdf(0, mean_life_length, 2.5));  % defines the survivability function of the cells
Z_normal=numerical_solution(S, f, 1-S, T, h);   % finds the solution for truncated normal life length G(t)

% calculates the solution for uniformly distributed life length
S=1-unifcdf(0:h:T, 0, mean_life_length*2)';  % defines the survivability function of the cells
Z_unif=numerical_solution(S, f, 1-S, T, h);

% calculates the solution for beta distributed life length
S=[1-betacdf(0:h/(2*mean_life_length):1, 2, 2)'; zeros(T/h-length(0:h/(2*mean_life_length):1)+1, 1)];
Z_beta=numerical_solution(S, f, 1-S, T, h);

% calculates the solution for chi^2 distributed life length
S=1-chi2cdf(0:h:T, mean_life_length)';
Z_chi2=numerical_solution(S, f, 1-S, T, h);

%% make plots of the results
line_wd=2.5;
% plot the different life length distributions used in the example
figure('visible','on', 'Units','pixels','OuterPosition',[0 0 1280 1024]);
set(gca,'FontSize',16)
hold on
plot(0:h:T, exppdf(0:h:T, 20)*h, 'Color', [1, 0, 0, 0.5], 'LineWidth', line_wd)
plot(0:h:T, normpdf(0:h:T, 10, 2.5)*h, '-', 'Color', [0, 0, 0, 0.5], 'LineWidth', line_wd)
plot(0:h:T, unifpdf(0:h:T, 0, 20)*h, '-', 'Color', [0, 1, 0, 0.5], 'LineWidth', line_wd)
plot(0:h:20, betapdf(0:(h/20):1, 2, 2)*h/20, '-', 'Color', [0, 0, 1, 0.5], 'LineWidth', line_wd)
plot(0:h:T, chi2pdf(0:h:T, 10)*h, '--', 'Color', [1, 0, 0, 0.5], 'LineWidth', line_wd)
xlabel('Time'); ylabel('Probability');
legend('Exp(10) life length', 'N(10, 2.5) life length', 'U(0, 20) life length', 'Beta(2,2) life length, scaled in [0, 20]', 'Chi^2(10) life length')
title('Life length distributions used in the example')
print('./figures/Example2_fig1', '-dpng', '-r0')

% plot the probability that successful mutant has not arrived yet for different life length distributions
figure('visible','on', 'Units','pixels','OuterPosition',[0 0 1280 1024]);
set(gca,'FontSize',16)
hold on
plot(0:h:T, Z_exp, 'Color', [1, 0, 0, 0.5], 'LineWidth', line_wd)
plot(0:h:T, Z_normal, '-', 'Color', [0, 0, 0, 0.5], 'LineWidth', line_wd)
plot(0:h:T, Z_unif, '-', 'Color', [0, 1, 0, 0.5], 'LineWidth', line_wd)
plot(0:h:T, Z_beta, '-', 'Color', [0, 0, 1, 0.5], 'LineWidth', line_wd)
plot(0:h:T, Z_chi2, '--', 'Color', [1, 0, 0, 0.5], 'LineWidth', line_wd)
xlabel('Time'); ylabel('Probability P(T>t)');
legend('Exp(10) life length', 'N(10, 2.5) life length', 'U(0, 20) life length', 'Beta(2,2) life length, scaled in [0, 20]', 'Chi^2(10) life length')
title('Successful mutant arrival in two-type Bellman-Harris BP')
print('./figures/Example2_fig2', '-dpng', '-r0')

% plot the distribution density function of T (time of successful mutant arrival)
figure('visible','on', 'Units','pixels','OuterPosition',[0 0 1280 1024]);
set(gca,'FontSize',16)
hold on
plot(h:h:T, -diff(Z_exp)*h, 'Color', [1, 0, 0, 0.5], 'LineWidth', line_wd)
plot(h:h:T, -diff(Z_normal)*h, '-', 'Color', [0, 0, 0, 0.5], 'LineWidth', line_wd)
plot(h:h:T, -diff(Z_unif)*h, '-', 'Color', [0, 1, 0, 0.5], 'LineWidth', line_wd)
plot(h:h:T, -diff(Z_beta)*h, '-', 'Color', [0, 0, 1, 0.5], 'LineWidth', line_wd)
plot(h:h:T, -diff(Z_chi2)*h, '--', 'Color', [1, 0, 0, 0.5], 'LineWidth', line_wd)
xlabel('Time'); ylabel('Probability Density Function of T');
legend('Exp(10) life length', 'N(10, 2.5) life length', 'U(0, 20) life length', 'Beta(2,2) life length, scaled in [0, 20]', 'Chi^2(10) life length')
title('P.D.F. of T for two-type Bellman-Harris BP')
text(40,4e-5,'Offspring p.g.f. for type 1 is f_1(s)=0.6250 + 0.375s^2', 'FontSize', 14)
print('./figures/Example2_fig3', '-dpng', '-r0')

