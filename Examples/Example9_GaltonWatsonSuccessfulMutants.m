%% Example 9: Calculating the distribution of successful mutant arrival in two-type Galton-Watson Branching Process (GWBP)
% This example pays respect to the first proposed model in which "successful mutants" are defined and presented.
% It considers the model in Serra, M.C., 2006. On the waiting time to escape. J. Appl. Probab. 43 (1),
% 296–302 and the theoretical continuations in Serra, M.C. and Haccou, P., 2007. Dynamics of escape mutants.
% Theoretical Population Biology, Vol 72:1, 167-178, https://doi.org/10.1016/j.tpb.2007.01.005
% The parameters of the model in this example are taken from the second article.
%
% Short description of the article model:
% The model considers Galton-Watson branching process - the life lengths are fixed, with length 1. 
% type 0 particles - supercritical, do not mutate to other types once they occur
% type 1 particles - subcritical, with probability of mutation to type 0
% The population starts from 1 particle of type 1.
% "Successful mutant" - a particle of type 0 that starts a type 0 branching, that never extincts
% "Unsuccessful mutant" - a particle of type 0 that starts a type 0 branching, that extincts
% The object of interest in the article is to calculate the distribution T of the first arrival of a "successful mutant".
%
% The distribution of the arrival time satisfies the same integral equation as in Examples 2 and 3, however substituting
% G(t)=0 for [0,1) and G(t)=1 for [1, +\infty) we get the results in the articles of Serra:
% P(T>k) = f_1(u*q_0+(1-u)*P(T>k-1)) dG(y), where the time k is discrete and integer. What this means is that we can
% actually use the integral equation and the numerical method for it with the only difference that G is a bit weird, 
% concentrated on a point. As the BHBP with exponential life length and the GWBP are both Markovian, the results are
% actually similar for them.
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
T=50;      % defines the interval [0, T] in which we want to find numerical solution

% defines the functions used in the integral equation
h=1;            % time is discrete with step 1
G=[0 1 ones(1, T-1)]';       % the life length distribution of Galton-Watson process, G must have length T/h+1=T+1
z=1-G;
q_0=0.5;
u=0.001;

%% binary splitting case
H=[1-0.375, 0, 0.375]';     % binary splitting that defines reproduction mean m=0.75, like in the article
f_1=@(x)(arrayfun(@(x)(sum(H'.*x.^(0:size(H,1)-1))), x));   % the offspring p.g.f. for type 1 cells
f=@(y, Z)(f_1(u*q_0+(1-u)*Z));  

P_successful_mutant=numerical_solution(z, f, G, T, h);

% we also have P(T>k, Z_1 (k)=0) = f(uq_0+(1-u)P(T>k, Z_1 (k)=0)), so:
P_joint=numerical_solution(zeros(size(G)), f, G, T, h);

% the hazard function, defined in the article is g(k)= P(T=k)/(P(T>k-1) - P(T>k-1, Z_1 (k-1)=0))
g=-diff(P_successful_mutant)./(P_successful_mutant(1:end-1)-P_joint(1:end-1));

%% Poisson distribution of the offspring
H=poisspdf(0:10, 0.75)';    % Poisson p.g.f. with lambda=0.75
f_1=@(x)(arrayfun(@(x)(sum(H'.*x.^(0:size(H,1)-1))), x));   % the offspring p.g.f. for type 1 cells
f=@(y, Z)(f_1(u*q_0+(1-u)*Z));  

P_mutant_pois=numerical_solution(z, f, G, T, h);

% we also have P(T>k, Z_1 (k)=0) = f(uq_0+(1-u)P(T>k, Z_1 (k)=0)), so:
P_joint_pois=numerical_solution(zeros(size(G)), f, G, T, h);

% the hazard function, defined in the article is g(k)= P(T=k)/(P(T>k-1) - P(T>k-1, Z_1 (k-1)=0))
g_pois=-diff(P_mutant_pois)./(P_mutant_pois(1:end-1)-P_joint_pois(1:end-1));

%% Linear fractional distribution of the offspring: f_1(s)=1-b/(1-c) + bs/(1-cs)
% custom case is geometric distribution with mean 0.75: f_1(s)=p/(1-(1-p)s), p = 4/7
f_1=@(s)(4/(7-3*s));
f=@(y, Z)(4./(7-3*(u*q_0+(1-u).*Z)));  

P_mutant_LF=numerical_solution(z, f, G, T, h);

% we also have P(T>k, Z_1 (k)=0) = f(uq_0+(1-u)P(T>k, Z_1 (k)=0)), so:
P_joint_LF=numerical_solution(zeros(size(G)), f, G, T, h);

% the hazard function, defined in the article is g(k)= P(T=k)/(P(T>k-1) - P(T>k-1, Z_1 (k-1)=0))
g_LF=-diff(P_mutant_LF)./(P_mutant_LF(1:end-1)-P_joint_LF(1:end-1));

%% make plots of the results
line_wd=2.5;
% plot the probability that successful mutant has not arrived yet
figure('visible','on', 'Units','pixels','OuterPosition',[0 0 1280 1024]);
set(gca,'FontSize',16)
hold on
plot(0:h:T, 1 - P_successful_mutant, 'x', 'Color', [1, 0, 0, 0.5], 'LineWidth', line_wd)
plot(0:h:T, 1 - P_mutant_pois, 'o', 'Color', [0, 0, 0.7, 0.2], 'LineWidth', line_wd)
plot(0:h:T, 1 - P_mutant_LF, ':', 'Color', [0, 0.7, 0, 0.7], 'LineWidth', line_wd)
xlabel('Time'); ylabel('Probability');
legend('Binary splitting', 'Poisson distribution', 'Linear fractional (geometric dist.)', 'Location', 'SouthEast')
title('Successful mutant arrival in two-type Galton-Watson BP, $P(T \leq k)$', 'interpreter', 'latex')
print('./figures/Example9_fig1', '-dpng', '-r0')

% The joint probability
figure('visible','on', 'Units','pixels','OuterPosition',[0 0 1280 1024]);
set(gca,'FontSize',16)
hold on
plot(0:h:T, P_joint, 'x', 'Color', [1, 0, 0, 0.5], 'LineWidth', line_wd)
plot(0:h:T, P_joint_pois, 'o', 'Color', [0, 0, 0.7, 0.5], 'LineWidth', line_wd)
plot(0:h:T, P_joint_LF, ':', 'Color', [0, 0.7, 0, 0.5], 'LineWidth', line_wd)
xlabel('Time'); ylabel('Probability');
legend('Binary splitting', 'Poisson distribution', 'Linear fractional (geometric dist.)', 'Location', 'SouthEast')
title('Successful mutant arrival in two-type Galton-Watson BP, $P(T > k, Z^1(k)=0)$', 'interpreter', 'latex')
print('./figures/Example9_fig2', '-dpng', '-r0')

% The hazard function g for binary splitting, poisson process and linear fractional (geometric distribution)
% note that in the article the legend is mixed up!
% so "poisson" is actually binary splitting
% "binary splitting" is actually linear fractional
% "linear fractional" is actually poisson
figure('visible','on', 'Units','pixels','OuterPosition',[0 0 1280 1024]);
set(gca,'FontSize',16)
hold on
plot(h:h:T, g, 'x', 'Color', [1, 0, 0, 0.5], 'LineWidth', line_wd)
plot(h:h:T, g_pois, 'o', 'Color', [0, 0, 0.7, 0.5], 'LineWidth', line_wd)
plot(h:h:T, g_LF, '^', 'Color', [0, 0.7, 0, 0.5], 'LineWidth', line_wd)
xlabel('Time'); ylabel('Hazard function g(k)');
legend('Binary splitting', 'Poisson distribution', 'Linear fractional (geometric dist.)', 'Location', 'SouthEast')
title('$g(k)= P(T=k)/[P(T>k-1) - P(T>k-1, Z_1 (k-1)=0)]$', 'interpreter', 'latex')
print('./figures/Example9_fig3', '-dpng', '-r0')


