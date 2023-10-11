function [xn,func_xn] = secantmethod(f,x0,x1,n)
% 
% Objective: Perform the secant algorithm method to approximate the root
% of a function with a given number of steps
%
% note: this function does not check whether the convergence criterion is
% met
%
% Input variables:
%   f - function handle whose root is desired
%   x0 - first of the two initial guesses for the root
%   x1 - second of the two initial guesses for the root
%   n - number of steps the secant method will perform
%
% Output variables:
%   xn - vector of the approximate roots of the function inputted at each
%       step
%   func_xn - values of the function evaluated at the approximate roots
%       (the values of xn)
%
%
% Functions called: none
%
%
% Create the xn vector and set the first two elements as x0 and x1
xn(1)=x0;
xn(2)=x1;
%
% Iterate through until xn will have (n+1) elements (since the first 
% element is x0)
for k=2:1:n
    %
    % Determine the next element of xn by finding where the secant line of
    % the previous two elements intersect with the x-axis
    xn(k+1)=xn(k)-(f(xn(k))*(xn(k)-xn(k-1)))/(f(xn(k))-f(xn(k-1)));
end
%
% Display the vector of the approximate roots of the function as well as
% the vector of the function evaluated at these approximate roots
xn
func_xn = f(xn)

end