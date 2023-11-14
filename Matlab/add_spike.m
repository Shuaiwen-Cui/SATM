function S_spike = add_spike(S,lb,ub)
%  S - signal
%  lb - lower boundary of spike magnitude
%  ub - upper boundary of spike magnitude
len = length(S);
time_spike = randi(len);
mag_spike = sign(rand(1)-0.5)*unifrnd(lb,ub);
S(time_spike) = mag_spike;
S_spike = S;