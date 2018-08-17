%% Example 4: Solving Renewal Equations and estimating the renewal functions
% Here we solve renewal equation presented in Example 2 from the Appendix of the article "Trayanov. P, Crump-Mode-Jagers
% Branching Process: A Numerical Approach; Eds: del Puerto, I. M., González, M., Gutiérrez, C., Martínez, R., Minuesa,
% C., Molina, M., Mota, M. and Ramos, A.; Lecture Notes in Statistics - Proceedings (219); Springer International 
% Publishing: 2016, pp. 167-182."
%
% The example solves the equation Z(t) = e^-t + \int_0^t Z(t-y) dG(y), for exponential distribution G(t) = 1-e^lambda*t.
% It is shown in the article that the theoretical solution is Z(t) = e^(-t) + lambda * (1 - e^(-t)).
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
T=10;           % defines the interval [0, T] in which we want to find numerical solution
h=0.001;        % defines the grid size of the numerical method

lambda=1/10;    % we set the average renewal time to be 10
G=expcdf(0:h:T, 1/lambda)';   % defines the distributions of the times between renewals to be exponential
z=exp(-(0:h:T));
f=@(y, Z)(Z);   % defines the equation as renewal type equation

% the expected number of renewals satisfies this equation
Z_numerical=numerical_solution(z, f, G, T, h);    % finds the solution, i.e. the expected population count in a Bellman-Harris branching process
Z_theoretical=exp(-(0:h:T)) + lambda * (1 - exp(-(0:h:T)));     % presented in the article

% make plots of the results
line_wd=2.5;
figure('visible','on', 'Units','pixels','OuterPosition',[0 0 1280 1024]);
set(gca,'FontSize',16)
hold on
plot(0:h:T, Z_numerical, 'Color', [0, 0, 0.7, 0.5], 'LineWidth', line_wd)
plot(0:h:T, Z_theoretical, '--', 'Color', [0, 0, 0, 0.5], 'LineWidth', line_wd)
xlabel('Time'); ylabel('Z(t)');
legend({'Numerical solution', 'Theoretical solution: $Z(t) = e^{-t} + \lambda (1 - e^{-t})$'}, 'interpreter', 'latex')
title('$Z(t) = e^{-t} + \int_0^t Z(t-y) dG(y)$, for exponential G with $\lambda=0.1$', 'interpreter', 'latex')
print('./figures/Example5_fig1', '-dpng', '-r0')




