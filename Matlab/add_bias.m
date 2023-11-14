function S_bias = add_bias(S,lb,ub)
%  S - signal
%  lb - lower boundary of bias magnitude
%  ub - upper boundary of bias magnitude
alb = abs(lb);
aub = abs(ub);
lb = min(alb,aub);
ub = max(alb,aub);
len = length(S);
time_bias = randi(len);
mag_bias = sign(rand(1)-0.5)*unifrnd(lb,ub);
for i = time_bias:len
    S(i) = S(i) + mag_bias;
end
S_bias = S;