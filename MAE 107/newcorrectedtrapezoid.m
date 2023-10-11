function I_CorrectedTRAP = newcorrectedtrapezoid(fx,a,b,n)

% Objective: Approximate the integral of a given function from a to b using
% corrected trapezoid rule with n steps
% 
% Input variables:
%   fx - function handle that will be subject to trapezoid rule
%   a - leftmost bound of the interval
%   b - rightmost bound of the interval
%   n - the number of steps used
%
% Output variables:
%   I_CTRAP - the computation of the approximate integral using corrected 
%       trapezoid rule
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
Trap = h*sum(f_xk)-(h/2)*(f_xk(1)+f_xk(n+1));

% Compute the approximate integral using corrected trapezoid by subtracting
% from the Trapezoid approximation according to the formula
I_CorrectedTRAP = Trap - (h/24)*(3*f_xk(1)-4*f_xk(2)+f_xk(3)+f_xk(n-1)-4*f_xk(n)+3*f_xk(n+1));

end
