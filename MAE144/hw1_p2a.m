% this defines a, b, and f using the RR_poly function
b = RR_poly([-2 2 -5 5], 1);
a = RR_poly([-1 1 -3 3 -6 6], 1);
f = RR_poly([-1 -1 -3 -3 -6 -6], 1);

% this runs the diophantine solver to find x and y such that a*x + b*y = f
[x, y] = RR_diophantine(a,b,f)

% this creates the expression for f using x,y from the diophantine solver
test_a = trim(a*x+b*y)
% this checks the difference between the expression for f using x,y from
%   the diophantine solver and the actual f that was initially defined 
residual1 = norm(f-test_a)

% Note that this answer gives an improper controller, as y is a 5th order
%    polynomial, while x is a 3rd order polynomial
% To fix this, we will change f to add additional poles to the controller
%   at s=-20 (see hw1_p2b)
