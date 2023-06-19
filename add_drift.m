function S_drift = add_drift(S,k_cap,SimStep)
%  S - signal
%  k_cap - max magnitude of drift slope
%  SimStep - time interval between frames
len = length(S);
k_cap = abs(k_cap);
time_drift = randi(len);
k = unifrnd(0,k_cap)*sign(rand(1)-0.5);
for i = time_drift:len
    S(i) = S(i) + k*(i-time_drift)*SimStep;
end
S_drift = S;