function [Dz_c, Dz_sc] = RPW_C2D_matched(bs,as,omegac,h)

% inputs: 'bs' is list of zeros, 'as' is list of poles, 'omegac' (optional)
%    is critical omega, 'h' is step size
% outputs: 'Dz_c' is the causal z-transform, 'Dz_sc' is the semi-causal
%   z-transform

% this is the mapping of zeros and poles from the s plane to the z plane
% for the matched z-transform method
map_zp_fun = @(s,h) exp(s*h);

% if only 3 input arguments, assume the third input is step-size
if nargin == 3
    h = omegac;
    omegac = 0;
end

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

% this creates the new zeros for semi-causal and causal systems
new_zeros_sc = -ones(1, number_of_inf_zeros);
new_zeros_c = -ones(1, number_of_inf_zeros-1);

% this adds the new zeros to the previously mapped zeros
bz_sc = cat(2, bz, new_zeros_sc);
bz_c = cat(2, bz, new_zeros_c);

% this formula calculates the gain of the matched z-transform for the
% critical omega of interest for both semi-causal and causal systems
k_sc = DCgain_star*(prod(-az + exp(1i*omegac*h))/prod(-bz_sc + exp(1i*omegac*h)));
k_c = DCgain_star*(prod(-az + exp(1i*omegac*h))/prod(-bz_c + exp(1i*omegac*h)));

% this makes z symbolic for future use
%sympref('FloatingPointOutput',true);
syms z;
% z = tf('z');

% this creates the numerator of the causal system by multiplying all of the
% (z-zeros) above
Dz_c_num = k_c*RR_Prod(z-bz_c);

% this creates the denominator of both systems by mutliplying all of the
% (z-poles) above
Dz_den = RR_Prod(z-az);

% this creates the causal D(z) z-transform
Dz_c = Dz_c_num/Dz_den

% this creates the numerator of the semi-causal system by multiplying all
% of the (z-zeros) 
Dz_sc_num = k_sc*RR_Prod(z-bz_sc);

% this creates the semi-causal D(z) z-transform
Dz_sc = Dz_sc_num/Dz_den

end