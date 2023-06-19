function Sf = Specdensity(ww,omegag,zetag,S0)
    Sf = S0*(omegag.^4+4*zetag.^2*omegag.^2.*ww.^2)./((omegag.^2-ww.^2).^2+4*zetag.^2*omegag.^2.*ww.^2);
end