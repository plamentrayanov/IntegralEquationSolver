%% Example 4: Solving Renewal Equations and estimating the renewal functions
% A renewal equation is integral equation of the type Z(t) = z(t) + \int_0^t Z(t-y) dG(y).
% In Renewal Theory its solution is represented by infinite sum of convolutions of increasing order:
% Z(t) = U*z(t), where U is the renewal function G^(0*)+ G^(1*) + ... + G^(n*) + ...
% This numerical method however does not calculate these convolutions numerically and is much faster.
% Also, the renewal function satisfies the equation U(t) = I(t) + \int_0^t U(t-y) dG(t) and we can find its solution 
% without calculating the convolutions. U(t)=N(t)+1, the number of renewals, counting the initial moment 0 as a renewal
% time.
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
T=100;          % defines the interval [0, T] in which we want to find numerical solution, i.e. the expected population count
h=0.01;         % defines the grid size of the numerical method

%% define a Poisson renewal process with averate rate 1/10
% the expected number of renewals (without counting the initial moment 0 as a renewal time) satisfies the equation
% N(t) = G(t) + \int_0^t N(t-y) dG(y).
% For Poisson process we know the solution is actually Z(t)=lambda*t, where lambda is the renewal rate (1/<the average 
% renewal interval length>). We can compare the theoretical and numerical solutions.
% The renewal function U(t) satisfies U(t) = I(t) + U * G

mean_life_length=10;    % we want to have an average renewal time of 10
G=expcdf(0:h:T, mean_life_length)';   % defines the distributions of the time between renewals to be exponential
z=G;
f=@(y, Z)(Z);   % defines the equation as renewal-type equation

% the expected number of renewals satisfies this equation
N_poisson_numerical=numerical_solution(z, f, G, T, h);    % estimates the solution, i.e. the expected number of renewals
N_poisson_theoretical=(0:h:T)/mean_life_length;

% the renewal function for Poisson process:
z=ones(T/h+1,1);
U_poisson_numerical=numerical_solution(z, f, G, T, h);  % estimates the renewal function as a solution to renewal equation
U_poisson_theoretical=1+(0:h:T)/mean_life_length;       % the theoretical solution is known in this example

% make plots of the results
line_wd=2.5;
figure('visible','on', 'Units','pixels','OuterPosition',[0 0 1280 1024]);
set(gca,'FontSize',16)
hold on
plot(0:h:T, N_poisson_numerical, 'Color', [1, 0, 0, 0.5], 'LineWidth', line_wd)
plot(0:h:T, N_poisson_theoretical, '--', 'Color', [0, 0, 0, 0.5], 'LineWidth', line_wd)
plot(0:h:T, U_poisson_numerical, 'Color', [0, 0, 0.7, 0.5], 'LineWidth', line_wd)
plot(0:h:T, U_poisson_theoretical, '--', 'Color', [0, 0, 0, 0.5], 'LineWidth', line_wd)
xlabel('Time'); ylabel('Expected number of renewals');
legend({'Numerical solution for $E(N(t))$', 'Theoretical Solution: $E(N(t))=\lambda t$', 'Numerical solution for $U(t)$', ...
    'Theoretical Solution: $U(t)=1+\lambda t$'}, 'interpreter', 'latex')
title('Poisson process. $N(t) = G(t) + N\ast G (t), U(t) = I(t) + U\ast G (t)$','interpreter', 'latex')
print('./figures/Example4_fig1', '-dpng', '-r0')

%% find the renewal function for uniform distribution G in [0, 1]
% in this case we know the theoretical solution for t in the interval [0, 1]: U(t) = e^t
T=1;
h=0.001;
G=unifcdf(0:h:T, 0, 1)';    % defines the distributions of the time between renewals to be exponential
f=@(y, Z)(Z);               % defines the equation as renewal type equation

z=ones(T/h+1,1);
U_uniform_numerical=numerical_solution(z, f, G, T, h);
U_uniform_theoretical=exp(0:h:T);

line_wd=2.5;
figure('visible','on', 'Units','pixels','OuterPosition',[0 0 1280 1024]);
set(gca,'FontSize',16)
hold on
plot(0:h:1, U_uniform_numerical, 'Color', [0, 0, 0.7, 0.5], 'LineWidth', line_wd)
plot(0:h:1, U_uniform_theoretical, '--', 'Color', [0, 0, 0, 0.5], 'LineWidth', line_wd)
xlabel('Time'); ylabel('Renewal Function');
legend({'Numerical solution for $U(t)$', 'Theoretical solution: $U(t)=e^t$ for $t \in [0,1]$'}, ...
    'Location', 'NW', 'interpreter', 'latex')
title('G is defined to be uniform in [0, 1]')
print('./figures/Example4_fig2', '-dpng', '-r0')



