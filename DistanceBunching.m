%{
==========================================================================
Bacis Information
FileName:       DistanceBunching.m
Author:         Shihao Zou
CreateTime:     2022.12.20
Version:        V2.0.0
Description:    This script is used for simulating distance bunching effect
in detection process of SWIM. Being constrained by the computing power,
normal laptop are not able to simulate this effect over large ocean area,
not even in the area scanned by SWIM (normally 18km X 18km). Thus the
simulated ocean size is limited to 6km X 6km, and grid size is 1m X 1m. The
setting of the rest parameters could be found in the codes. You can
minimize the value of each parameter, but be careful while increasing the
variables, it could cause overload of computing devices.

==========================================================================
%}

clear;
close all;

% Loading parameters from json config file.
Conf = jsondecode(fileread("ParaConfig.json"));
ExpCounts = Conf.ExpParameters.ExpCounts;
IncidentAngleIndex = Conf.ExpParameters.IncidentAngleIndex;
OceanWaveType = Conf.ExpParameters.OceanWaveType;
OceanParameters = Conf.ExpParameters.OceanParameters;
ObservingIntervalAngle = Conf.ExpParameters.ObservingIntervalAngle * pi/180;

IncidentAngle = Conf.RadarParameters.IncidentAngles(IncidentAngleIndex) * pi/180;
RadarGates = Conf.RadarParameters.RadarGates(IncidentAngleIndex);
RpetitionPulses = Conf.RadarParameters.PulsesEachMicroCycle(IncidentAngleIndex);
TransmitPower = Conf.RadarParameters.TransmitPower;
LightSpeed = Conf.RadarParameters.LightSpeed;
WorkingFrequency = Conf.RadarParameters.WorkingFrequency;
WorkingWavelength = LightSpeed/WorkingFrequency;
RotatingSpeed = Conf.RadarParameters.RotatingSpeed * pi/180;
RadarAltitude = Conf.RadarParameters.RadarAltitude;
% note that the size of the simulated ocean area would be equal to the
% project of Three3dBBeamWidth on the horizontal plane.

ObserveCountsPerCycle = floor(2 * pi / ObservingIntervalAngle);
% Observing internal time
delta_t = ObservingIntervalAngle / (RotatingSpeed / 60);

% ocean size.
Ly = 6000;
Lx = Ly;
Three3dBBeamWidth = atan(Ly * cos(IncidentAngle) / RadarAltitude);
% Parameter in Gaussian distrbution to make sure antenna gain at the 
% boundary of ocean area is half of the peak value.
Ly3dB = Ly / (2 * sqrt(2 * log(2)));

dx = 1;
dy = dx;
M = Lx / dx;
N = M;
if mod(M, 4) ~= 0
    M = 4 * floor(M / 4);
    N = M;
end
if M > 6000
    disp('Warning: Ocean spliting size over 6000 X 6000.');
end

% refresh ocean area size
Lx = M * dx;
Ly = N * dy;

disp(['Simulated ocean size is: ',num2str(Lx), 'mX', num2str(Ly), ...
    'm, Scanning interval is: ', num2str(ObservingIntervalAngle*180/pi), ...
    'o, Exp ', num2str(ExpCounts),  ' rounds'])

% ====================Simulation Process======================== %
OmniPmod = zeros(ObserveCountsPerCycle, M / 2);
OmniSlopeSpectrum = zeros(ObserveCountsPerCycle, M / 2);
OmniPmodNoised = zeros(ObserveCountsPerCycle, M / 2);
OmniSlopeSpectrumNoised = zeros(ObserveCountsPerCycle, M / 2);

for ObserveCount = 1 : ObserveCountsPerCycle
    OmniPmodNoised_tmp = 0;
    OmniSlopeSpectrumNoised_tmp = 0;
    OmniPmod_tmp = 0;
    OmniSlopeSpectrum_tmp = 0;
    for ExpRound=1 : ExpCounts

        %% Forward Process of SWIM detection, output receiving power.
        ObservingAngle = (ObserveCount - 1) * ObservingIntervalAngle;
        t = (ExpRound - 1) * delta_t;
        disp(['Start processing for Exp round: ',num2str(ExpRound), '. Observing angle: ', num2str(ObservingAngle*180/pi)]);

        % generating sea surface elevation, slope along x, and slope along
        % y axis.
        [height, slope_x, slope_y] = GeneratingSeaSurface(Lx, Ly, dx, dy, t, ObservingAngle);

        % calculating back coef averaged by anttenna gain pattern.
        [Sigma0, ~] = CalculatingAveragedBackCoef(height, slope_x, slope_y, Lx, Ly, dx, dy, ObservingAngle, Conf);
        
        % simulating the distance bunching effect process.
        [Sigma0] = RearrangeSigma(height, Lx, Ly, dx, dy, Sigma0, Conf);
        
        % first zero-centering the sigma and normalize it.
        [Sigma0_Relative, Sigma0_mean] = ZeroCenteredAndNormalization(Lx, Ly, dx, dy, RadarAltitude, IncidentAngle, Sigma0);

        % simulate speckle noise for each transmitted pulses. Note this
        % simulation method for speckle noise origins from Hauser et al.
        % 2001.
        ReceivingPower = 0;
        for c_pulse = 1 : RpetitionPulses
            [Sigma0_Relative_noised] = AddingSpeckleNoise(Sigma0_Relative + 1, RadarGates);

            [ReceivingPower0] = GeneratingReceivingPower(Lx, Ly, dx, dy, Sigma0_Relative_noised, Sigma0_mean, Conf);
            ReceivingPower = ReceivingPower + ReceivingPower0;
            disp(['*** The ',num2str(c_pulse), 'th pulses simulated.']);

        end
        ReceivingPower = ReceivingPower / RpetitionPulses;

        %% Inversion Process of SWIM detection, output receiving power.
        % calculate modulated signal from echo power
        [modulated_signal, alpha] = RemovingDeterminsticParameters(Lx, Ly, dx, dy, ReceivingPower, Conf);

        % inversing the modulation spectrum from modulated signal with
        % alpha factor
        % processing modulated signal with speckle noise.
        InversionMode = 0;
        [Pmod, SlopeSpectrum] = EstimateModulationSpectrum(Lx, Ly, dx, dy, modulated_signal, alpha, RpetitionPulses, IncidentAngle, InversionMode);
        OmniPmodNoised_tmp = OmniPmodNoised_tmp + Pmod;
        OmniSlopeSpectrumNoised_tmp = OmniSlopeSpectrumNoised_tmp + SlopeSpectrum;

        % processing modulated signal without speckle noise.
        InversionMode = 1;
        [Pmod, SlopeSpectrum] = EstimateModulationSpectrum(Lx, Ly, dx, dy, Sigma0_Relative, alpha, RpetitionPulses, IncidentAngle, InversionMode);
        OmniPmod_tmp = OmniPmod_tmp + Pmod;
        OmniSlopeSpectrum_tmp = OmniSlopeSpectrum_tmp + SlopeSpectrum;
    end
    
    % Exp round averaged Omni Modulation spectrum from noise-free and
    % noised modulated signal
    OmniPmod(ObserveCount, :) = OmniPmod_tmp / ExpCounts;
    OmniSlopeSpectrum(ObserveCount, :) = OmniSlopeSpectrum_tmp / ExpCounts;

    OmniPmodNoised(ObserveCount, :) = OmniPmodNoised_tmp / ExpCounts;
    OmniSlopeSpectrumNoised(ObserveCount, :) = OmniSlopeSpectrumNoised_tmp / ExpCounts;
end

K = linspace(0.0001, pi / dx, M / 2);
Theta = (0 : ObserveCountsPerCycle - 1) * ObservingIntervalAngle * 180 / pi;

figure;
imagesc(K, Theta, OmniSlopeSpectrum);
title("Inversed Slope Spectrum from noise-free signal in the Polar System.");
xlabel("Wavenumber k (1/m)");
ylabel("Azimuth (o)");

figure;
imagesc(K, Theta, OmniSlopeSpectrumNoised);
title("Inversed Slope Spectrum from noised signal in the Polar System.");
xlabel("Wavenumber k (1/m)");
ylabel("Azimuth (o)");

SlopeSpectrumAlongFlight = OmniSlopeSpectrum(1, :);
SlopeSpectrumNoisedAlongFlight = OmniSlopeSpectrumNoised(1, :);
figure;
plot(K, SlopeSpectrumAlongFlight);
hold on;
plot(K, SlopeSpectrumNoisedAlongFlight);
title("Comparison between slope Spectrum along flight directionl.");
xlabel("Wavenumber k (1/m)");
ylabel("SlopeSpectrum (m^2)");
legend("SlopeSpec-NoiseFree","SlopeSpec-Noised");

