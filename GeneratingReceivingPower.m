%{
==========================================================================
Bacis Information
FileName:       GeneratingReceivingPower.m
Author:         Shihao Zou
CreateTime:     2018.04
Version:        V1.0.0
Description:    This function is used calculate the receiving power of
radar at given sea conditions.

==========================================================================
Input parameters
Lx: sea surface size along x axis.
Ly: sea surface size along y axis.
dx: the size of spliting bin of sea surface along x axis.
dy: the size of spliting bin of sea surface along y axis.
Sigma0_Relative_noised: Noised value of the relativity of backscattering coef.
Sigma0_mean: Mean value of the averaged backscattering coef.
Conf: Configuration file content

Output parameters:
ReceivingPower: Receiving power of radar at given sea condition.
%}

function [ ReceivingPower ] = GeneratingReceivingPower(Lx, Ly, dx, dy, Sigma0_Relative_noised, Sigma0_mean, Conf)
    RadarAltitude = Conf.RadarParameters.RadarAltitude;
    IncidentAngleIndex = Conf.ExpParameters.IncidentAngleIndex;
    IncidentAngle = Conf.RadarParameters.IncidentAngles(IncidentAngleIndex)*pi/180;
    Three3dBBeamWidth = Conf.RadarParameters.Three3dBBeamWidth*pi/180;
    TransmitPower = Conf.RadarParameters.TransmitPower;
    LightSpeed = Conf.RadarParameters.LightSpeed;
    WorkingFrequency = Conf.RadarParameters.WorkingFrequency;
    WorkingWavelength = LightSpeed/WorkingFrequency;

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

    PositionVector=sqrt((x+RadarAltitude*tan(IncidentAngle)).^2+(y).^2+(RadarAltitude).^2);
    IncidentAnglesPerBin = RadarAltitude./PositionVector;
    IncidentAnglesCosPerBin=sum(IncidentAnglesPerBin')/N;
    R=RadarAltitude./IncidentAnglesCosPerBin';
    thetaR=acos(IncidentAnglesCosPerBin');
    betaf=Three3dBBeamWidth/2/sqrt(2*log(2));
    Ge_2=(exp(-((thetaR-IncidentAngle)./(sqrt(2)*betaf)).^2)).^2.*(26000/Three3dBBeamWidth)^2;%天线地距向增益的平方
    CR=(TransmitPower*WorkingWavelength./(4*pi.*R).^3).*Ge_2*dx;%雷达功率系数，请参考雷达回波功率方程公式

    Ly3dB = Ly/(2*sqrt(2*log(2)));
    Gy = exp(-y.^2/(2*Ly3dB^2));
    ReceivingPower=CR'.*Sigma0_mean.*Sigma0_Relative_noised*dy*sum(Gy(1,:).^2);
end