%% Example 7: Estimating the probability generating function (P.G.F.) of Bellman-Harris branching process
% This example is actually a special case of the previous Example7, as the p.g.f. of such process satisfies the 
% equation F(t,s) = s(1-G(t)) + \int_0^t f(F(t-u, s)) dG(u), where f(s) = \Sum_0^\infty p_n*s^n. 
% The probabilities p_n(u), for having n children, are constant, independent of the age.
% The offspring p.g.f. f(s) defines a supercritical branching process.
% TIP: ADJUST THE DEFINITION OF f(s) TO SEE RESULTS FOR SUBCRITICAL OR CRITICAL PROCESS!
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
T=50;               % defines the interval [0, T] in which we want to find numerical solution
delta=0.1;          % defines the grid size of the numerical method
s_delta=0.01;       % add fine refinement for s as we are interested in better approximation for s partial derivatives

G=1-(1-normcdf(0:delta:T, 10, 2.5)')./(1-normcdf(0, 10, 2.5));    % mean life length for women is 76
G(end)=1;   % death probability for ages greater than omega is 1

H=[0.3 0.4 0.2 0.1]';   % contains p_n - probabilities for having n children
f=@(u, s)(polyval(flipud(H), s));   % f depends only on s, in contrast to Example7, where f depends on u and s

% estimate F(t, s) for s in [0, 1]:
F=zeros(size(G));
s=0:s_delta:1;
parfor i=1:length(s)
    z=s(i)*(1-G);
    F(:,i)=numerical_solution(z, f, G, T, delta);
end

%% calculate the expected population as dF/ds(t, 1) with s_delta=0.01:
first_derivative=diff(F,1,2)./s_delta;
second_derivative=diff(F,2,2)./s_delta^2;
ExpPop=first_derivative(:,end);
VAR_Pop = second_derivative(:, end) + first_derivative(:,end) - first_derivative(:,end).^2;
VAR_Pop(VAR_Pop<0)=0;   % due to rounding error the calculated variation could be slightly below 0, so adjust it before taking the root
STD_Pop=sqrt(VAR_Pop);

%% make plots of the results
line_wd=2.5;
custom_colororder=[1:-1/(size(F,2)-1):0; zeros(2, size(F,2))]';

figure('visible','on', 'Units','pixels','OuterPosition',[0 0 1280 1024], 'defaultAxesColorOrder', custom_colororder);
set(gca,'FontSize',16)
hold on
h_plot=plot(0:delta:T, F, 'LineWidth', line_wd);
xlabel('Time'); ylabel('F(t,s) for s=0:0.01:1');
legend([h_plot(1), h_plot(end)], 'F(t,0) = P(Z(t) = 0) \rightarrow q_0', 'F(t,1) = 1 theoretically')
title('P.G.F. of Bellman-Harris BP with N(10, 2.5) life length')
text(10, 0.1, ['The probabilities for giving birth are independent of age for Bellman-Harris BP. \newline'...
    'The offsring p.g.f. is 0.3 + 0.4s + 0.2s^2 + 0.1s^3 for all ages, i.e. custom case of Sevastyanov BP.'],...
    'FontSize',14)
print('./figures/Example8_fig1', '-dpng', '-r0')

% make 3D plot of F(t,s) for Sevastyanov Process
figure('visible','on', 'Units','pixels','OuterPosition',[0 0 1280 1024]);
set(gca,'FontSize',16)
colormap(jet)
hold on
mesh(0:delta:T, 0:s_delta:1, F', 'FaceAlpha',0.5, 'FaceColor', 'interp');
material('SHINY')
lighting('GOURAUD')
light
view(3)
xlabel('t'); ylabel('s'); zlabel('F(t,s)')
title('P.G.F. of Bellman-Harris BP with N(10, 2.5) life length')
print('./figures/Example8_fig2', '-dpng', '-r0')

% Make 3D plot of dF/ds(t,s) for Sevastyanov Process
figure('visible','on', 'Units','pixels','OuterPosition',[0 0 1280 1024]);
set(gca,'FontSize',16)
colormap(jet)
hold on
mesh(0:delta:T, s_delta:s_delta:1, diff(F,1,2)', 'FaceAlpha',0.5, 'FaceColor', 'interp');
material('SHINY')
lighting('GOURAUD')
light
view(3)
xlabel('t'); ylabel('s'); zlabel('$\frac{sF}{ds}(t,s)$','interpreter','latex')
title('Bellman-Harris BP with N(10, 2.5) life length. Partial derivative on s: $\frac{dF}{ds}(t,s)$','interpreter','latex')
text(1, 0.9, 0.02, ['The probabilities for giving birth are independent of age for Bellman-Harris BP. \newline'...
    'The offsring p.g.f. is 0.3 + 0.4s + 0.2s^2 + 0.1s^3 for all ages, i.e. custom case of Sevastyanov BP.'],...
    'FontSize', 14)
print('./figures/Example8_fig3', '-dpng', '-r0')

% Draw the expected population as dF/ds(t,1)
figure('visible','on', 'Units','pixels','OuterPosition',[0 0 1280 1024], 'defaultAxesColorOrder', custom_colororder);
set(gca,'FontSize', 16)
hold on
plot(0:delta:T, ExpPop, 'LineWidth', line_wd);
plot(0:delta:T, STD_Pop, '--', 'LineWidth', line_wd);
xlabel('Time'); ylabel('Expectation and Volatility');
legend({'$E(Z)=\frac{dF}{ds}(t,1)$', '$\sigma(Z)=\sqrt{\frac{d^2F}{ds^2}(t,1) + \frac{dF}{ds}(t,1)-\frac{dF^2}{ds}(t,1)}$'},...
    'interpreter', 'latex')
text(20, 0.5, 'The offsring p.g.f. is 0.3 + 0.4s + 0.2s^2 + 0.1s^3.', 'FontSize', 14)
title('Central moments of Bellman-Harris BP with N(10, 2.5) life length')
print('./figures/Example8_fig4', '-dpng', '-r0')

