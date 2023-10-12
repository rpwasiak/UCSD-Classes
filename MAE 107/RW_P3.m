%
% Objective: this script uses 4th order Runge-Kutta method with varying 
%   number of steps to solve the problem y' = g(t,y) with initial value 
%   y(0) = 4 where g(t,y) is the solution of 
%   z = exp(-(1+sin(z)))-(sin(t+y)^2)*((1+z^2)^(1/3)) which is solved using
%   fixed point method
%
%
%
% Functions called: 
%   fixedpointeval - this function 
% 
%
%
% Set final T value to 6
T = 6;
%
% Set function handle of g(t,y,z) (which is set = z)
g = @(t,y,z) exp(-(1+sin(z)))-(sin(t+y)^2)*((1+z^2)^(1/3));
%
% Iterate through 7 times
for i=1:7
    % Set number of steps, n, to 4*2^i
    n = 2^(i+2);
    % Create n_vector to store number of steps used for each iteration
    n_vector(i) = n;
    % Set step-size, h, equal to T/n
    h = T/n;
    %
    % Set first index of each row of t equal to 0
    t(i,1) = 0; 
    % Set first index of each row of yRK4 equal to 4
    yRK4(i,1) = 4;
    %
    % Iterate through n times
    for j=1:n
        %
        % RK-4th order
        % Determine K1,K2,K3, and K4 values using the 4th order Runge-Kutta 
        %   formula where f(t,y) is the solution of z = g(t,y,z), which is 
        %   determined using fixedpointeval with 10 steps where t and y are
        %   the respective index in their matrices, and initial z_guess=1
        K1RK4 = h*fixedpointeval(g,1,t(i,j),yRK4(i,j),10);
        K2RK4 = h*fixedpointeval(g,1,t(i,j)+(h/2),yRK4(i,j)+K1RK4/2,10);
        K3RK4 = h*fixedpointeval(g,1,t(i,j)+(h/2),yRK4(i,j)+K2RK4/2,10);
        K4RK4 = h*fixedpointeval(g,1,t(i,j)+h,yRK4(i,j)+K3RK4,10);
        %
        % Evaluate the next value of the row using the RK4 formula
        yRK4(i,j+1) = yRK4(i,j)+(1/6)*(K1RK4+2*(K2RK4+K3RK4)+K4RK4);
        %
        % Set the next value in the i'th row of t equal to h + current t value
        t(i,j+1) = t(i,j)+h;
        %
    end
    %
    % Create a plot on figure 1 of the respective yRK4 values vs the i'th 
    %   row of t (up to the n+1 index) with a title, axis labels, and a 
    %   legend (for the different number of steps used)
    figure(1);
    plot(t(i,[1:n+1]),yRK4(i,[1:n+1]));
    title("4th order Runge Kutta Method of solution of y = g(t,y)");
    title(legend(num2str(n_vector'),'Location','southwest'),'Number of Steps');
    ylabel('y');
    xlabel('t');
    hold on
    %
end
%
% This function is the fixed-point method for a function g that takes 3
%   inputs (t, y, and z) that is set equal to z and outputs z_sol which is
%   the n+1 index of z
function z_sol = fixedpointeval(g,z0,t0,y0,n)
% Set first index of z to z0
z(1) = z0;
%
% Iterate through n times
for k=1:n
    % Set next value of z equal to g evaluated at t0, y0, and current z
    %   value
    z(k+1) = g(t0,y0,z(k));
end
% Set z_sol equal to the last index of z
z_sol = z(end);
%
end