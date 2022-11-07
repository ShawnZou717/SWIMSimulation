%{
==========================================================================
Bacis Information
FileName:       ZeroCenteredAndNormalization.m
Author:         Shihao Zou
CreateTime:     2018.04
Version:        V1.0.0
Description:    This function is zero-centering the averaged back coef and
then normalize it.

==========================================================================
Input parameters
Lx: sea surface size along x axis.
Ly: sea surface size along y axis.
dx: the size of spliting bin of sea surface along x axis.
dy: the size of spliting bin of sea surface along y axis.
RadarAltitude
IncidentAngle
Sigma0: Averaged back coef by anttenna gain pattern.

Output parameters:
Sigma0_Relative: Backscattering coefficients Relativity
Sigma0_mean: Mean value of back coef Sigma0
%}

function [ Sigma0_Relative, Sigma0_mean ] = ZeroCenteredAndNormalization( Lx, Ly, dx, dy, RadarAltitude, IncidentAngle, Sigma0 )
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
    
    x=tan(sqrt(1-IncidentAnglesCosPerBin.^2)./IncidentAnglesCosPerBin).^2;
    y=log(Sigma0)+log(IncidentAnglesCosPerBin.^4);
    P=polyfit(x,y,1);
    Sigma0_mean=polyval(P,x);
    Sigma0_mean = exp(Sigma0_mean-log(IncidentAnglesCosPerBin.^4));
    Sigma0_Relative=(Sigma0-Sigma0_mean)./Sigma0_mean;
    Sigma0_Relative=Sigma0_Relative;
    Sigma0_mean=Sigma0_mean;
end