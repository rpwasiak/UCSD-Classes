function I_LE = newleftendpoint(fx,a,b,n)

% Objective: Approximate the integral of a given function from a to b using
% left-endpoint rule with n steps
% 
% Input variables:
%   fx - function handle that will be subject to left-endpoint rule
%   a - leftmost bound of the interval
%   b - rightmost bound of the interval
%   n - the number of steps used
%
% Output variables:
%   I_LE - the computation of the approximate integral using left-endpoint 
%       rule
%
% Functions called;
%   evalf - evaluates f(x) at x in increments of (b-a)/n from a to b
%
%

% Compute the step-size, h
h = (b-a)/n;

% Use the evalf function to create a vector of f(x) evaluated at x in
% increments of h from a to b
f_xk = evalf(fx,a,b,n);


% Compute the approximate integral using left-endpoint rule (note: we
% subtract f_xk(n+1) because this value is f(b) which is not used in
% left-endpoint rule)
I_LE = h*(sum(f_xk) - f_xk(n+1));

end