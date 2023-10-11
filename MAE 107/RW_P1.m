%
% Objective: this script uses varying number of steps while employing
%   Euler's method, two Second Order Runge-Kutta methods (with B2 = 1/2 and
%   with B2 = 1/3), and Fourth Order Runge-Kutta method to determine the 
%   the solution up to T=2 of y' = e^(-y)*sin(t+2*pi*y) with initial 
%   condition y(0)=1 
%
%
%
% Functions called: none
%   
% 
%
% Set function handle of f(t,y) = y'
f = @(t,y) exp(-y)*sin(t+2*pi*y);
% Set final value T = 2
T = 2;
%
% Find the approximate solution of y(t) using RK4 for 8192 steps
% Set the first index of est_solution = 1
est_solution(1) = 1;
% Set the first index of t_star = 0
t_star(1) = 0;
% Set the step size, h_star, as T/8192
h_star = T/8192;
%
% Iterate through 8192 times (for 8192 steps)
for a = 1:8192
    % Set the K1,K2,K3, and K4 values using the respective formulas for the
    %   RK4 method
    K_1RK4 = h_star*f(t_star(a),est_solution(a));
    K_2RK4 = h_star*f(t_star(a)+h_star/2,est_solution(a)+K_1RK4/2);
    K_3RK4 = h_star*f(t_star(a)+h_star/2,est_solution(a)+K_2RK4/2);
    K_4RK4 = h_star*f(t_star(a)+h_star,est_solution(a)+K_3RK4);
    % 
    % Set the next value of est_solution using the formula with K1,K2,K3
    %   and K4
    est_solution(a+1) = est_solution(a)+(1/6)*(K_1RK4+2*(K_2RK4+K_3RK4)+K_4RK4);
    % 
    % Set the next value of t_star equal to h + current t_star value
    t_star(a+1) = t_star(a)+h_star;
end
%
% Iterate through 11 times
for i=1:11
    % Set n equal to 2*2^i
    n = 2^(i+1);
    % Create n_vector to store number of steps used for each iteration
    n_vector(i) = n;
    % Set step-size, h, equal to T/n
    h = T/n;
    %
    % Set first index of each row of t equal to 0
    t(i,1) = 0; 
    % Set first index of each row of yEuler, yRK2_1, yRK2_2, and yRK4 equal
    %   to 1
    yEuler(i,1) = 1;
    yRK2_1(i,1) = 1;
    yRK4(i,1) = 1;
    yRK2_2(i,1) = 1;
    %
    % Iterate through n times
    for j=1:n
        % 
        % Perform Euler's: Set next index of i'th row equal to current
        %   index of i'th row + stepsize*f(t,y)
        yEuler(i,j+1) = yEuler(i,j)+h*f(t(i,j),yEuler(i,j));
        %
        %
        %
        % Perform RK-2nd order with B2 = 1/2
        % Set B2_half = 1/2 and use the three equations to find B1_half,
        %   VRK2 and ORK2
        B2_half = 1/2;
        B1_half = 1 - B2_half;
        VRK2 = 1/(2*B2_half);
        ORK2 = 1/(2*B2_half);
        %
        % Set the K1 and K2 values using the formulas from Second Order RK 
        %   method with h, f, VRK2, and ORK2
        K1RK2 = h*f(t(i,j),yRK2_1(i,j));
        K2RK2 = h*f(t(i,j)+(h*VRK2),yRK2_1(i,j)+(ORK2*K1RK2));
        %
        % Set the next index of the i'th row of yRK2_1 using the Second
        %   Order RK formula with B1, B2, K1, and K2
        yRK2_1(i,j+1) = yRK2_1(i,j)+B1_half*K1RK2+B2_half*K2RK2;
        %
        %
        %
        % Peform RK-2nd order with B2 = 1/3
        % Set B2_new = 1/2 and use the three equations to find B1_new,
        %   VRK2_2 and ORK2_2
        B2_new = 1/3;
        B1_new = 1 - B2_new;
        VRK2_2 = 1/(2*B2_new);
        ORK2_2 = 1/(2*B2_new);
        %
        % Set the K1 and K2 values using the formulas from Second Order RK 
        %   method with h, f, VRK2_2, and ORK2_2
        K1RK2_2 = h*f(t(i,j),yRK2_2(i,j));
        K2RK2_2 = h*f(t(i,j)+(h*VRK2_2),yRK2_2(i,j)+(ORK2_2*K1RK2_2));
        %
        % Set the next index of the i'th row of yRK2_2 using the Second
        %   Order RK formula with B1_new, B2_new, K1, and K2
        yRK2_2(i,j+1) = yRK2_2(i,j)+B1_new*K1RK2_2+B2_new*K2RK2_2;
        %
        %
        %
        % Perform RK-4th order
        % Set the K1,K2,K3, and K4 values using the respective formulas for the
        %   RK4 method 
        K1RK4 = h*f(t(i,j),yRK4(i,j));
        K2RK4 = h*f(t(i,j)+(h/2),yRK4(i,j)+K1RK4/2);
        K3RK4 = h*f(t(i,j)+(h/2),yRK4(i,j)+K2RK4/2);
        K4RK4 = h*f(t(i,j)+h,yRK4(i,j)+K3RK4);
        %
        % Set the next value of yRK4 using the formula with K1,K2,K3 and K4
        yRK4(i,j+1) = yRK4(i,j)+(1/6)*(K1RK4+2*(K2RK4+K3RK4)+K4RK4);
        %
        %
        %
        % Set the next value in the i'th row of t equal to h + current t value
        t(i,j+1) = t(i,j)+h;
        %
    end
    %
    % Create a plot on figure 1 of the respective yEuler values vs the i'th 
    %   row of t (up to the n+1 index) with a title, axis labels, and a 
    %   legend (for the different number of steps used)
    figure(1);
    plot(t(i,[1:n+1]),yEuler(i,[1:n+1]));
    title("Euler's Method of e^-^y*sin(t+2*pi*y)");
    title(legend(num2str(n_vector'),'Location','northwest'),'Number of Steps');
    ylabel('y');
    xlabel('t');
    hold on
    %
    % Create a plot on figure 2 of the respective yRK2_1 values vs the i'th 
    %   row of t (up to the n+1 index) with a title, axis labels, and a 
    %   legend (for the different number of steps used)
    figure(2);
    plot(t(i,[1:n+1]),yRK2_1(i,[1:n+1]));
    title("2nd order Runge Kutta Method with B2 = 1/2 of e^-^y*sin(t+2*pi*y)");
    title(legend(num2str(n_vector'),'Location','northwest'),'Number of Steps');
    ylabel('y');
    xlabel('t');
    hold on
    %
    % Create a plot on figure 3 of the respective yRK2_2 values vs the i'th 
    %   row of t (up to the n+1 index) with a title, axis labels, and a 
    %   legend (for the different number of steps used)
    figure(3);
    plot(t(i,[1:n+1]),yRK2_2(i,[1:n+1]));
    title("2nd order Runge Kutta Method with B2 = 1/3 of e^-^y*sin(t+2piy)");
    title(legend(num2str(n_vector'),'Location','northwest'),'Number of Steps');
    ylabel('y');
    xlabel('t');
    hold on
    %
    % Create a plot on figure 4 of the respective yRK4 values vs the i'th 
    %   row of t (up to the n+1 index) with a title, axis labels, and a 
    %   legend (for the different number of steps used)
    figure(4);
    plot(t(i,[1:n+1]),yRK4(i,[1:n+1]));
    title("4th order Runge Kutta Method of e^-^y*sin(t+2piy)");
    title(legend(num2str(n_vector'),'Location','northwest'),'Number of Steps');
    ylabel('y');
    xlabel('t');
    hold on
    %
    %
    % Create the error vectors for Eulers, RK2_1, RK2_2, and RK4 and set
    %   the i'th index as the absolute value of the difference between the
    %   value y(T), which is est_solution(end), and the n+1 index of the
    %   i'th row of each method
    error_euler(i) = abs(yEuler(i,n+1)-est_solution(end));
    error_RK2_1(i) = abs(yRK2_1(i,n+1)-est_solution(end));
    error_RK2_2(i) = abs(yRK2_2(i,n+1)-est_solution(end));
    error_RK4(i) = abs(yRK4(i,n+1)-est_solution(end));
    %
end
%
% Create 4 plots on figure 5 of the log(error vectors) vs 
%   log(respective number of steps used) with a title, axis labels, and a
%   legend (for the different steps used)
figure(5);
title('Log-Log plot of Errors vs Number of steps');
xlabel('Log10(n)');
ylabel('Log10(error of y(T))');
hold on
plot(log10(n_vector),log10(error_euler));
plot(log10(n_vector),log10(error_RK2_1));
plot(log10(n_vector),log10(error_RK2_2));
plot(log10(n_vector),log10(error_RK4));
title(legend({"Euler's",'RK2 with B2 = 1/2','RK2 with B2 = 1/3','RK4'}),'Method Used');
%
%
%
%
%
% Problem 2
%
% Create 4 plots on figure 6 of the log(error vectors) vs 
%   log(respective number of function evaluations used) with a title, axis 
%   labels, and a legend (for the different steps used)
figure(6);
title('Log-Log plot of Errors vs Number of Function Evaluations');
xlabel('Log10(F)');
ylabel('Log10(error of y(T))');
hold on
% Euler's method uses 1 function evaluation per step, so n steps = n
%    function evaluations
plot(log10(n_vector),log10(error_euler));
% RK2 method uses 2 function evaluation per step, so n steps = 2n
%    function evaluations (hence the 2*n_vector)
plot(log10(2*n_vector),log10(error_RK2_1));
plot(log10(2*n_vector),log10(error_RK2_2));
% RK4 method uses 4 function evaluation per step, so n steps = 4n
%    function evaluations (hence the 4*n_vector)
plot(log10(4*n_vector),log10(error_RK4));
title(legend({"Euler's",'RK2 with B2 = 1/2','RK2 with B2 = 1/3','RK4'}),'Method Used');



