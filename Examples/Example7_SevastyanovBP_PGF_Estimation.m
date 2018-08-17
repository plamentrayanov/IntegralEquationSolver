%% Example 7: Estimating the probability generating function (P.G.F.) of Sevastyanov branching process
% The p.g.f. of such process satisfies the equation F(t,s) = s(1-G(t)) + \int_0^t f(u, F(t-u, s)) dG(u), where
% f(u, s) = \Sum_0^\infty p_n(u)*s^n. The probabilities p_n(u), for having n children, depend on the age u.
% The offspring p.g.f. f(u,s) is defined at the end of the script!
%
% TIP: The example presents a supercritical branching process. To see the results for subcritical or critical branching
% process, adjust the definition of f(u,s).
%
% From the p.g.f F(t,s) we can find its derivatives and calculate the Expected population count and its volatility.
% There is also an alternative way for calculating the expectated population count of the branching process:
% The p.g.f. of such process satisfies the equation F(t,s) = s(1-G(t)) + \int_0^t f(u, F(t-u, s)) dG(u).
% After differentiating we find an integral equation for the expected population count:
% E(Z(t))=dF/ds(t,1) = 1-G(t) + \int_0^t (df/du(u,1) + df/dv(u,1))EZ(t-u) dG(u).
% Alternative approach is, instead of estimating F(t,s) and dF/ds(t,1), we can directly solve the equation for E(Z(t)).
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
T=50;           % defines the interval [0, T] in which we want to find numerical solution
delta=0.1;      % defines the grid size of the numerical method
s_delta=0.01;   % add fine refinement for s as we are interested in better approximation for s partial derivatives

G=1-(1-normcdf(0:delta:T, 10, 2.5)')./(1-normcdf(0, 10, 2.5));    % mean life length is 10
G(end)=1;   % death probability for ages greater than omega is 1

% estimate F(t, s) for s in [0, 1]:
F=zeros(size(G));
s=0:s_delta:1;
parfor i=1:length(s)
    z=s(i)*(1-G);
    F(:,i)=numerical_solution(z, @f, G, T, delta);
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

% Make 2D plot of F(t,s) for Sevastyanov Process
figure('visible','on', 'Units','pixels','OuterPosition',[0 0 1280 1024], 'defaultAxesColorOrder', custom_colororder);
set(gca,'FontSize',16)
hold on
h_plot=plot(0:delta:T, F, 'LineWidth', line_wd);
xlabel('Time'); ylabel('F(t,s) for s=0:0.01:1');
legend([h_plot(1), h_plot(end)], 'F(t,0) = P(Z(t) = 0) \rightarrow q_0', 'F(t,1) = 1 theoretically')
title('P.G.F. of Sevastyanov BP with N(10, 2.5) life length')
text(15, 0.1, ['the probabilities for giving birth change with age \newline'...
    'the offsring p.g.f. changes from 0.7 + 0.1s + 0.1s^2 + 0.1s^3 at age 0\newline', ...
    'to 0.3 + 0.4s + 0.2s^2 + 0.1s^3 at age 10 and stays constant for ages greater than 10'],...
    'FontSize',14)
print('./figures/Example7_fig1', '-dpng', '-r0')

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
title('P.G.F. of Sevastyanov BP with N(10, 2.5) life length')
print('./figures/Example7_fig2', '-dpng', '-r0')

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
title('Sevastyanov BP with N(10, 2.5) life length. Partial derivative on s: $\frac{dF}{ds}(t,s)$','interpreter','latex')
text(0, 0.9, 0.014, ['The probabilities for giving birth depend on the age for Sevastyanov BP. \newline'...
    'The offsring p.g.f. changes from 0.7 + 0.1s + 0.1s^2 + 0.1s^3 at age 0\newline', ...
    'to 0.3 + 0.4s + 0.2s^2 + 0.1s^3 at age 10 and stays constant for ages greater than 10.'],...
    'FontSize',14)
print('./figures/Example7_fig3', '-dpng', '-r0')

% Draw the expected population as dF/ds(t,1)
figure('visible','on', 'Units','pixels','OuterPosition',[0 0 1280 1024], 'defaultAxesColorOrder', custom_colororder);
set(gca,'FontSize',16)
hold on
plot(0:delta:T, ExpPop, 'LineWidth', line_wd);
plot(0:delta:T, STD_Pop, '--', 'LineWidth', line_wd);
xlabel('Time'); ylabel('Expectation and Volatility');
legend({'$E(Z)=\frac{dF}{ds}(t,1)$', '$\sigma(Z)=\sqrt{\frac{d^2F}{ds^2}(t,1) + \frac{dF}{ds}(t,1)-\frac{dF^2}{ds}(t,1)}$'},...
    'interpreter', 'latex')
title('Central moments of Sevastyanov BP with N(10, 2.5) life length')
print('./figures/Example7_fig4', '-dpng', '-r0')

%% h definition
function res=f(u, s)
% defines the offspring p.g.f. as a function of the age u.
% the probabilities for giving birth change with age - the offsring p.g.f. move from 0.7 + 0.1s + 0.1s^2 + 0.1s^3 at 
% age 0, to 0.3 + 0.4s + 0.2s^2 + 0.1s^3 at age 10 and stay like this for age greater than 10
% The mean number of children for the first law is 0.6, which defines a subcritical process, but for the second law is
% 1.1, which defines a supercritical process. The law depends on the age of the individual and it moves linearly from
% the subcritical process for young ages to the supercritical for older ages. The mean life length is 10.
res=zeros(length(u),1);
for i=1:length(res)
    H=(min(10,u(i))*[0.7 0.1 0.1 0.1]' + max(10-u(i), 0)*[0.3 0.4 0.2 0.1]')/10;
    res(i)=polyval(flipud(H), s(i));
end
end

