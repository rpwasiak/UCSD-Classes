function [a0,a1,a2] = least_squares_quadratic(points)
%
% Objective: this function performs the least squares method for a 
%   quadratic fit: x(t) = a0 + a1*x + a2*x^2
%
% 
% Input variables: 
%   points - this is a matrix of 2 columns with rows being (t,x) of as many
%       points as the user has
%
% Output variable:
%   a0 - the coefficient for the constant of x(t)
%   a1 - the coefficient for x of x(t)
%   a2 - the coefficient for x^2 of x(t)
%
%
% Functions called: none
%  
%
% Determine the number of points the user inputted
r = size(points,1);
%
% Iterate through for each point
for k=1:r
    %
    % Take the t and x values from the input matrix and put them into the
    %   vector of t and x
    t(k)=points(k,1);
    x(k)=points(k,2);
    % 
    % Evaluate the t^4, t^3, t^2, (t^2)*x, and t*x values for the current
    %   point, which is used for the least squares method for quadratic
    %   approximation (comes from minimizing difference between x(t_i) and 
    %   x_i))
    c1(k) = t(k)^4;
    c2(k) = t(k)^3;
    c3(k) = t(k)^2;
    c4(k) = (t(k)^2)*x(k);
    c5(k) = t(k)*x(k);
end
%
% Create the A and b matrix used in the least squares method by summing the
%   values computed above and placing them in the right spot of the matrix
A = [sum(c1) sum(c2) sum(c3); sum(c2) sum(c3) sum(t); sum(c3) sum(t) r]
b = [sum(c4); sum(c5); sum(x)]
%
% Solve the Ax=b matrix problem, which yields x_t = [a2;a1;a0]
x_t = A\b;
%
% Assign the values from x_t to its corresponding value (a0,a1,or a2)
a0 = x_t(3)
a1 = x_t(2)
a2 = x_t(1)
%
% Display the least squares quadratic approximation, x(t) to the user
xt = ['x(t) = ',num2str(a0),' + ',num2str(a1),'x + ',num2str(a2),'x^2'];
disp(xt)
%
% Plot the Data Points and the x(t) approximation, title the plot, label
% the axes and create a legend
hold on
i = 0:0.01:8;
f = a0+a1*i+a2*i.^2;
plot(i,f)
plot(t,x,'r*')
xlabel('t')
ylabel('x(t)')
title('Least Squares Method for Quadratic Estimate')
legend('Least Squares Quadratic Approximation','Data Points')
%
%
end