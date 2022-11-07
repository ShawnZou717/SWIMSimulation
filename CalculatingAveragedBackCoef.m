%{
==========================================================================
Bacis Information
FileName:       CalculatingBackCoef.m
Author:         Shihao Zou
CreateTime:     2018.04
Version:        V1.0.0
Description:    This function is used to generate Backscattering 
coefficients of the generated sea surface.

==========================================================================
Input parameters
height: sea surface elevations.
slope_x: slope of sea surface along x axis.
slope_y: slope of sea surface along y axis.
Lx: sea surface size along x axis.
Ly: sea surface size along y axis.
dx: the size of spliting bin of sea surface along x axis.
dy: the size of spliting bin of sea surface along y axis.
ObservingAngle
Conf: Configuration file content

Output parameters:
Sigma0: Backscattering coefficients of each bin of sea surface.
%}

function [ Sigma0 ] = CalculatingAveragedBackCoef( height, slope_x, slope_y, Lx, Ly, dx, dy, ObservingAngle, Conf)
RadarAltitude = Conf.RadarParameters.RadarAltitude;
IncidentAngleIndex = Conf.ExpParameters.IncidentAngleIndex;
IncidentAngle = Conf.RadarParameters.IncidentAngles(IncidentAngleIndex)*pi/180;
OceanWaveType = Conf.ExpParameters.OceanWaveType;
OceanParameters = Conf.ExpParameters.OceanParameters;

% the number of spliting bins of sea surface
M=floor(Lx/dx);
N=floor(Ly/dy);

% index matrix in wavenumber domain
m=0:M-1;
m=m'*ones(1,M);
n=0:N-1;
n=ones(1,N)'*n;

% index matrix in Cartesian space domain
p=m;
x = -Lx/2+p*dx;
q=n;
y = -Ly/2+q*dx;

% calculating the angle between the position vector from surface point
% to spectrometer and the normal vector of ocean section on each ocean
% bin.
PositionVector=abs(slope_x.*(x+RadarAltitude*tan(IncidentAngle))+slope_y.*y+1*(RadarAltitude-height));
MagnitudePositionVector=sqrt((x+RadarAltitude*tan(IncidentAngle)).^2+(y).^2+(RadarAltitude-height).^2);
MagnitudeNormalVector=sqrt((slope_x).^2+(slope_y).^2+(1).^2);
AngleCosValue = PositionVector./(MagnitudePositionVector.*MagnitudeNormalVector);
AngleTanValue = sqrt(1-AngleCosValue.^2)./AngleCosValue;

% Backscattering Coef is calculated based on Cox-Munk coefficient model
lo=0.61;

WindDirection = nan;
for index=1:length(OceanWaveType)
    if strcmp(OceanWaveType(index), "EL") || strcmp(OceanWaveType(index), "PM")...
        || strcmp(OceanWaveType(index), "JSON")
         ocean_paras = OceanParameters(index);
         ocean_paras = ocean_paras{1,1};
         WindDirection = ocean_paras.PeakWaveDirection*pi/180;
         WindSpeed = ocean_paras.WindSpeed;
    end
end

if isnan(WindDirection)
    error('There is no wind wave in the ocean type.')
end

cs_fine_local=cos(ObservingAngle-WindDirection);
sn_fine_local=sqrt(1-cs_fine_local.^2);

% Empirical euqation between slope variance and win speed, given by the
% doctoral dissertation of Xiaoqing Chu.
su=sqrt(0.00078545*WindSpeed+0.0092407);
sc=sqrt(0.00052799*WindSpeed+0.0097295);

v=2./((cs_fine_local./su).^2+(sn_fine_local./sc).^2);
pp=exp(-AngleTanValue.^2./v)./(2*pi*su*sc);

Sigma0=lo*pi*AngleCosValue.^(-4).*pp;

% sigma0 averaged with attenna pattern
Gy = exp(-y.^2/(2*Ly^2));
Gy2=(Gy').^2;
Sigma0=Gy2.*Sigma0';
Sigma0=sum(Sigma0)./sum(Gy2);

end