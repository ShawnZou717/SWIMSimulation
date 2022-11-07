function [ Pmod, SlopeSpectrum ] = EstimateModulationSpectrum(Lx, Ly, dx, dy, modulated_signal, alpha, RpetitionPulses, IncidentAngle, InversionMode)
    M = Lx/dx;
    
    mfft = fft(modulated_signal);
    mfft=fftshift(mfft);
    Pm = mfft.*conj(mfft)/M^2*Lx/2/pi;
    Pm = Pm(end-M/2+1:end);

    k=linspace(0.0001,pi/dx,M/2);
    kp=2*sqrt(2*log(2))/(0.47/sin(IncidentAngle));
    Psp=1/(sqrt(2*pi)*kp*RpetitionPulses).*exp(-(k.*k)/(2*kp*kp));
    
    if InversionMode==0
        Pmod = Pm-Psp;
        SlopeSpectrum = Pmod/(sqrt(2*pi)/Ly*mean(alpha)^2);
    elseif InversionMode==1
        Pmod = Pm;
        SlopeSpectrum = Pmod/(sqrt(2*pi)/Ly*mean(alpha)^2);
    end
end