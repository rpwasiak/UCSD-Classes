function xn = fixedpointmethod(g,x0,n)
% Objective: this function performs the fixed-point method to determine the
%   solution of the form g(x) = x
% 
% Input variables: 
%   g - function handle that is set equal to x (ie. we want to solve
%       g(x) = x)
%   x0 - initial guess of the solution
%   n - arbitrary number of steps picked by the user
%
% Output variable:
%   fsolutions - the vector of the approximate solution, x, at each step
%
%
% Functions called: none
% 
%
%
% Create the xn vector and set the first element as the initial guess, x0
xn(1) = x0;
%
% Iterate through n times
for k=1:1:n
    % 
    % Determine the next element of xn by evaluating g(x) at the value of 
    % the current xn element
    xn(k+1) = g(xn(k));
end

end