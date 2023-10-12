function I_TRAP = newtrapezoid(fx,a,b,n)

% Objective: Approximate the integral of a given function from a to b using
% trapezoid rule with n steps
% 
% Input variables:
%   fx - function handle that will be subject to trapezoid rule
%   a - leftmost bound of the interval
%   b - rightmost bound of the interval
%   n - the number of steps used
%
% Output variables:
%   I_TRAP - the computation of the approximate integral using trapezoid
%       rule
%
% Functions called:
%   evalf - evaluates f(x) at x in increments of (b-a)/n from a to b
%
%

% Compute the step-size, h
h = (b-a)/n;

% Use the evalf function to create a vector of f(x) evaluated at x in
% increments of h from a to b
f_xk = evalf(fx,a,b,n);

% Compute the approximate integral using trapezoid rule 
% note: we subtract (h/2)*(f_xk(1)+f_xk(n+1) because this value is 
% (h/2)(f(a)+f(b)) and we multiplied h(f(a)+f(b)) so subtracting this makes
% sure we are only adding (h/2)(f(a)+f(b))
I_TRAP = h*sum(f_xk)-(h/2)*(f_xk(1)+f_xk(n+1));

end