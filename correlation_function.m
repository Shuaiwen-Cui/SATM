function R = correlation_function(t,S0,omige_g,zeta_g)
    %金井清功率谱密度函数
    b1 = omige_g*(1+4*zeta_g^2)/zeta_g;
    b2 = omige_g*(1-4*zeta_g^2)/sqrt(1-zeta_g^2);
    omige_d = omige_g*sqrt(1-zeta_g^2);
    R = 0.5*pi*S0*exp(-zeta_g*omige_g*abs(t)).*(b1*cos(omige_d*t)+b2*sin(omige_d*abs(t)));
end