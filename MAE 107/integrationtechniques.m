function [realeval,ILE,ITRAP,ICTRAP] = integrationtechniques(fun,a,b,k) 

% Objective: Use left-endpoint, trapezoid, and corrected trapezoid
% rules to approximate the integral of a given function using different 
% numbers of steps
% 
% Also, find the errors of the integration techniques with varying numbers
% of steps and plot these errors on the same graph
%
% Input variables:
%   fun - function handle that will be subject to integration techniques
%   a - left bound of the interval
%   b - right bound of the interval
%   k - 10, 10^2, 10^3 ... 10^k steps will be used
%
% Output variables:
%   realeval - the true solution
%   ILE - vector of the integral approximations using left-endpoint with
%       increasing number of steps
%   ITRAP - vector of the integral approximations using trapezoid with
%       increasing number of steps
%   ILE - vector of the integral approximations using corrected trapezoid
%       with increasing number of steps
%
%   A log-log plot of the errors of each technique vs the number of steps
%   will also be produced
%
%
% Functions called:
%   newleftendpoint - computes the integral of the function using
%       left-endpoint rule
%   newtrapezoid - computes the integral of the function using trapezoid
%       rule
%   newcorrectedtrapezoid - computes the integral of the function using
%       corrected trapezoid rule
%   (all of these functions also use evalf - yields a vector of f(x) values
%   at increments of h over the interval from a to b)


% Compute the true solution of the integral using the "integral" function
realeval = integral(fun,a,b)

% Iterate through and create vectors of the integral approximations using
% left-endpoint, trapezoid, and corrected trapezoid rules with 10^i steps
for i=1:k
    ILE(i) = newleftendpoint(fun,a,b,10^i);
    ITRAP(i) = newtrapezoid(fun,a,b,10^i);
    ICTRAP(i) = newcorrectedtrapezoid(fun,a,b,10^i);
    step(i) = 10^i;
end

% Display the vectors of the integral approximations of left-endpoint,
% trapezoid, and corrected trapezoid found using the varying steps
ILE
ITRAP
ICTRAP

% Create vectors of the errors of each integration approximation for the 
% different techniques at the different step sizes
eLE = abs(ILE - realeval);
eTRAP = abs(ITRAP - realeval);
eCTRAP = abs(ICTRAP - realeval);

% Plot the results on a log-log plot, with the x-axis: log10(steps used)
% and the y-axis: log10(error of each integral approximation)
% Plot is also properly titled, labeled, and includes a legend
figure(1);
hold on
plot(log10(step),log10(eLE),'LineWidth',4)
plot(log10(step),log10(eTRAP),'LineWidth',2)
plot(log10(step),log10(eCTRAP),'LineWidth',2)
xlabel('Log10(n)');
ylabel('Log10(errors)');
title('Errors of Integral Approximations of f(x)');
legend('Error of Left-endpoint rule','Error of Trapezoid rule','Error of Corrected Trapezoid rule');
end
