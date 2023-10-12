function f_xk = evalf(fx,a,b,n)

% Objective: Evaluate a given function f(x) at x in increments of (b-a)/n 
% over the interval from a to b and output this vector
% 
% Input variables:
%   fx - function handle that will be evaluated at different x values
%   a - leftmost bound of the interval
%   b - rightmost bound of the interval
%   n - the number of increments desired
%
% Output variables:
%   f_xk - a vector of f(x) evaluated at x in increments of (b-a)/n over
%       the interval from a to b
%   
%

% Compute the size of the increment, h
h = (b-a)/n;

% Create a vector from 0 to n
k = 0:1:n;

% Create a vector from a to b, with increments of h
x_k = (a + h.*k);

% Create a vector of f(x) evaluated at values of x over the interval from a
% to b, with increments of h
f_xk = fx(x_k);

end