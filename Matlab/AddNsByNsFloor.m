function NsSig = AddNsByNsFloor(SigClean,NsLevel)
% SigClean: Signal to be added noise
% NsLevel: Unit dB
% NoiseBase: usually a standard Gaussian White Noise
seg = length(SigClean);
NoiseBase = idinput(seg,'RGS',[0,1],[-1,1]); % standard noise as reference for noise
rms_ns_bs = rms(NoiseBase);
amplify_fac = sqrt(exp(NsLevel/10));
rms_ns = amplify_fac * rms_ns_bs;
Ns = rms_ns * NoiseBase;
NsSig = SigClean + Ns;
end