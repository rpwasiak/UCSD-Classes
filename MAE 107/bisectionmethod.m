function fzero = bisectionmethod(g,a0,b0,eps)
% 
% Objective: Perform the bisection algorithm to approximate the root of a
% function inside an interval within a desired error
%
% note: this function does not check whether there is indeed a root inside
% the interval (doesn't check for bad input data)
%
% Input variables:
%   g - function handle whose root is desired
%   a0 - left bound of interval
%   b0 - right bound of interval
%   eps - desired max error of root approximation
%
% Output variable:
%   fzero - approximate root of the function inputted
%
%
% Functions called: none
%
%
%
% Determine the number of steps based on the size of the interval and the
% max error desired and rounding the step count up to the nearest integer 
n=ceil(log2((b0-a0)/eps));
%
% Create the b vector and initialize the first element as b0
b(1)=b0;
% Create the a vector and initialize the first element as a0
a(1)=a0;
% Create the c vector and initialize the first element as the midpoint 
% between a0 and b0
c(1)=(a0+b0)/2;
%
% Iterate n times the following code
for k=1:1:n
    % Determine the  next element of c as the midpoint between the last a
    % and b elements
    c(k+1)=(a(k)+b(k))/2;
    %
    % Check if the function evaluated at the midpoint of new interval
    % equals 0 and if so, set fzero equal to the midpoint, and end the 
    % for loop
    if g(c(k+1))==0
        fzero = c(k+1);
        break
    % Multiply the function evaluated at this new midpoint with the
    % function evaluated at the current right bound and if this is less
    % than 0, make the next left bound this midpoint and keep the right 
    % bound the same. If this is not less than zero, make the next right
    % bound this midpoint and keep the left bound the same
    elseif g(c(k+1))*g(b(k)) < 0
        a(k+1)=c(k+1);
        b(k+1)=b(k);
    else
        b(k+1)=c(k+1);
        a(k+1)=a(k);
    end
    % Check if the loop is on its last step, and if it is, set fzero equal
    % to the current new midpoint
    if k==n
        fzero = c(n+1);
    end
end

end




