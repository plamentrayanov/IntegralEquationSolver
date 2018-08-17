function Z=numerical_solution(z, f, G, T, h)
%       numerical_solution - opensource program for solving integral equations of the type 
%       Z(t) = z(t) + \int_0^t [f(y, Z(t-y))] dG(y)
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
%
% Z=numerical_solution(z, f, G, T, h) solves numerically the equation Z(t)=z(t) + \int_0^t f(y, Z(t-y))dG(y) and 
% returns a vector Z with the estimated values on the grid 0:h:T.
%
% z and G are vectors with length ceil(T/h)+1, which represent the values of the functions at points 0:h:T.
%
% f is a function of two arguments.
%
% T defines the interval [0, T] in which we estimate the solution.
%
% h determines the grid size of the numerical method. The smaller it is, the better the estimation error.
%
% Note: z, f and G are input functions that determine the equation uniquely. The method works for smooth functions with 
% finite number of jump type discontinuities.

Z=zeros(ceil(T./h)+1,1);
Z(1)=z(1);
for k=2:length(Z)
    % it is important to use the transposed (1:k-1)', as it allows the use of two-parameter arrayfun to define f
    Z(k)=z(k)+ f((1:k-1)', Z(k-1:-1:1))'*(G(2:k)-G(1:k-1));
end
end