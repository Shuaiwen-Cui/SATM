function [EqData,EqTime,TailVal] = EqSimulation(SimStep,EqDur,omegag,zetag,S0,fac_time,pos_time1,pos_time2)
% simulated by using Kanai-Tajimi model and superposing harmonic waves
% [1] SimStep (sec) - time interval
% [2] T (sec) - duration of earthquake
% [3] omegag - parameter in KT model
% [4] zetag - parameter in KT model
% [5] S0 - base of specdensity curve in KT model
% [6] fac_time - factor controlling the time modulating curve
% [7] pos_time1 - parameter to control the time modulating curve
% [8] pos_time2 - parameter to control the time modulating curve


N0 = 200; % 200-600
sample = 1; % one earthquake sample each time
a = -pi; % lower limit of random variables
b = pi; % upper limit of random variables
dt = SimStep;
EqDur = floor(EqDur/dt)*dt; % Duration 
t = dt:dt:EqDur;
wu = 1/SimStep/2;   %upper limit of frequency
dw = wu/N0;
ww = [dw:dw:wu]; % frequency vector

gt = fac_time * (exp(-pos_time1*t)-exp(-pos_time2*t)); % time modulating

Ag1 = zeros(length(t),length(ww));

for i=1:sample
    theda1=a+(b-a)*randn(1,sample);
end

for j = 1:sample
    for k = 1:N0
        Xk0(k) = sqrt(2)*cos(k*theda1(j)+pi/4);
        Yk0(k) = sqrt(2)*sin(k*theda1(j)+pi/4);
    end
    Xk(j,:) = Xk0(randperm(N0,N0));
    Yk(j,:) = Yk0(randperm(N0,N0));
end

k = [1:1:N0];
for j=1:sample
    for i=1:length(t)
        AA=sqrt(2*SpecDensity(dw.*k,omegag,zetag,S0)*dw).*(cos(dw.*k*i*dt).*Xk(j,k)+sin(dw.*k*i*dt).*Yk(j,k));
        Ag2(i)=sum(AA);
    end
    Ag(j,:) = Ag2.*gt;
end

EqTime = t;
EqData = Ag(1,:);
TailVal = abs(gt(end));
% figure
% plot(EqTime,gt,'r-',EqTime,EqData,'g-')
end



