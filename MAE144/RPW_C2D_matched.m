function [Dz] = RPW_C2D_matched(bs,as,omegac,h,causality)

% inputs: 'bs' is list of zeros, 'as' is list of poles, 'omegac' (optional)
%    is critical omega, 'h' is step size, 'causality' is 1 for causal and 0
%    for semi-causal
% outputs: 'Dz' is the z-transform, with desired causality

% this is the mapping of zeros and poles from the s plane to the z plane
% for the matched z-transform method
map_zp_fun = @(s,h) exp(s*h);

% this performs the mapping if there are zeros
if ~isempty(bs)
    bz = map_zp_fun(bs,h);
else 
    bz = bs;
end

% this performs the mapping if there are poles
if ~isempty(as)
    az = map_zp_fun(as,h);
else
    az = as;
end

% this calculates the gain of D(s) at the critical omega of interest
DCgain_star = prod(-bs+(1i*omegac))/prod(-as+(1i*omegac));

% this determines the number of infinite zeros being added according to the
% matched z-transform method
number_of_inf_zeros = length(as) - length(bs);

% this creates the new zeros for desired semi-causal or causal system
new_zeros = -ones(1, number_of_inf_zeros-causality);

% this adds the new zeros to the previously mapped zeros
bz_new = cat(2, bz, new_zeros);

% this formula calculates the gain of the matched z-transform for the
% critical omega of interest
k = DCgain_star*(prod(-az + exp(1i*omegac*h))/prod(-bz_new + exp(1i*omegac*h)));

% this makes z symbolic for future use
%sympref('FloatingPointOutput',true);
syms z;
% z = tf('z');

% this creates the numerator of the z-transform by multiplying all of the
% (z-zeros) above
Dz_num = k*RR_Prod(z-bz_new);

% this creates the denominator of the z-transform by mutliplying all of the
% (z-poles) above
Dz_den = RR_Prod(z-az);

% this creates the D(z) z-transform
Dz = Dz_num/Dz_den


end