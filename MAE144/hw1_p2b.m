% this defines a, b, and f_new with 6 additional poles at s=-20 using the RR_poly function
b = RR_poly([-2 2 -5 5], 1);
a = RR_poly([-1 1 -3 3 -6 6], 1);
f_new = RR_poly([-1 -1 -3 -3 -6 -6 -20 -20 -20 -20 -20 -20], 1);

% this runs the diophantine solver to find x and y s.t a*x + b*y = f_new
[x_new, y_new] = RR_diophantine(a,b,f_new)

% these two lines check the difference between the expression for f using x,y from
%   the diophantine solver and the actual f_new that was initially defined 
test_b = trim(a*x_new+b*y_new)
residual2 = norm(f_new-test_b)

% Note that now, the numerator, y, is still a 5th order polynomial, but the
% denominator, x, is now a 6th order polynomial, making this controller
% proper

