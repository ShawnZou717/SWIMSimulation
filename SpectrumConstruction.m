%{
==========================================================================
Bacis Information
FileName:       SpectrumConstruction.m
Author:         Shihao Zou
CreateTime:     2018.04
Version:        V1.0.0
Description:    This function is used to generate sea surface with given
sea wave model, such as PL model, DV model and EL model.

==========================================================================
Input parameters
Lx: sea surface size along x axis.
Ly: sea surface size along y axis.
dx: the size of spliting bin of sea surface along x axis.
dy: the size of spliting bin of sea surface along y axis.
ObservingAngle: the angle between the scanning direction and flight
direction.
ocean_type:  Type of wave model such as "DV", "PM"
ocean_paras: specific parameters for generating given type of ocean waves

Output parameters:
F0: wave spectrum of given ocean type.
%}

function [FD] = SpectrumConstruction(Lx, Ly, dx, dy, ObservingAngle, ocean_type, ocean_paras)
% the number of spliting bins of sea surface
M=floor(Lx/dx);
N=floor(Ly/dy);

% resolution in wavenumber domain
dkx=2*pi/Lx;
dky=2*pi/Ly;

% index matrix in wavenumber domain
m=0:M;
m=m'*ones(1,M+1);
n=0:N;
n=ones(1,N+1)'*n;
kx = -pi/dx+m*dkx;
ky = -pi/dy+n*dky;
K = sqrt(kx.^2+ky.^2);

% acceleration speed of gravity
g0=9.81;

if strcmp(ocean_type,'PM')
    WindSpeed = ocean_paras.WindSpeed;
    WindDirection = ocean_type.PeakWaveDirection*pi/180;
    k_peak=0.7694*g0/(WindSpeed);

    delta_phi=zeros(M+1,N+1);
    delta_phi(1:M/2,N/2+1:N+1)=pi;
    delta_phi(1:M/2,1:N/2)=pi;
    delta_phi(M/2+1:M+1,1:N/2)=2*pi;

    fine_transmit=delta_phi+atan((ky./kx));
    FD=exp(-5/4*(k_peak./K).^2)*(0.0081/2).*K.^(-4).*4.*cos((WindDirection-ObservingAngle)-fine_transmit).^4./(3*pi);%PM谱波高谱变量

elseif strcmp(ocean_type,'JONS')
    WindSpeed = ocean_paras.WindSpeed;
    WindDirection = ocean_type.PeakWaveDirection*pi/180;
    kafang=90000;%60km风区
    kafang_w=g0*kafang/WindSpeed;
    alfa=0.076*(kafang_w)^(-0.22);
    gama=3.3;
    w0_w=18.3*kafang_w^(-1/3);
    w0=w0_w*g0/WindSpeed;
    k_peak=w0^2/g0;

    deltas_l1=zeros(M+1,N+1);
    deltas_s1=zeros(M+1,N+1);
    [row1,col1]=find(K>k_peak);
    for hh1=1:length(row1)
        deltas_l1(row1(hh1),col1(hh1))=0.09;
    end
    [row2,col2]=find(K<=k_peak);
    for hh1=1:length(row2)
        deltas_s1(row2(hh1),col2(hh1))=0.07;
    end
    deltas1=deltas_l1+deltas_s1;
    delta_phi=zeros(M+1,N+1);
    delta_phi(1:M/2,N/2+1:N+1)=pi;
    delta_phi(1:M/2,1:N/2)=pi;
    delta_phi(M/2+1:M+1,1:N/2)=2*pi;

    fine_transmit=delta_phi+atan((ky./kx));
    FD=alfa/2.*K.^(-4).*exp(-5/4.*(k_peak./K).^2).*gama.^exp(-((K./k_peak).^0.5-1).^2./(2.*deltas1.^2)).*4.*cos((WindDirection-ObservingAngle)-fine_transmit).^4./(3*pi);

elseif strcmp(ocean_type,'DV')
    WaveHeight = ocean_paras.WaveHeight;
    PeakWavelength = ocean_paras.PeakWavelength;
    PeakWaveDirection = ocean_paras.PeakWaveDirection*pi/180;
    deltas=0.006;
    Hs=WaveHeight;
    k_peak=2*pi/PeakWavelength;
    FD_=Hs.^2./(16*sqrt(2*pi)*deltas).*exp(-0.5.*((K-k_peak)/deltas).^2);
    
    
    syms xx yy;
    cp=cos(xx-yy).^14;
    cpp=int(cp,xx,0,2*pi);
    
    F322D1=cos((PeakWaveDirection-ObservingAngle)).^14./double(cpp);

    delta_phi=zeros(M+1,N+1);
    delta_phi(1:M/2,N/2+1:N+1)=pi;
    delta_phi(1:M/2,1:N/2)=pi;
    delta_phi(M/2+1:M+1,1:N/2)=2*pi;

    fine_transmit=delta_phi+atan((ky./kx));
    F322D=cos((PeakWaveDirection-ObservingAngle)-fine_transmit).^14./double(cpp);
    FD=FD_./K.*F322D;

elseif strcmp(ocean_type,'EL')
    WindSpeed = ocean_paras.WindSpeed;
    InverseWaveAge = ocean_paras.InverseWaveAge;
    WindDirection = ocean_paras.PeakWaveDirection*pi/180;
    omega = InverseWaveAge;
    km = 363;
    cm = sqrt(2*g0/km); 
    ap = 0.006.*sqrt( omega );
    k0 = g0 ./ (WindSpeed.^2);
    kp = k0 .* omega.^2;
    cp = WindSpeed/omega;
    sigma = 0.08*(1+4*omega.^(-3));
    [us]=blpara(WindSpeed);
        
    if us<cm
        am = 0.01*(1+log(us/cm));
    elseif us>cm
        am = 0.01*(1+3*log(us/cm));
    end
    if omega>=0.84 && omega<1
        gama = 1.7;
    elseif omega>1 && omega<5
        gama = 1.7 + 6*log10(omega);
    end

    cc1 = sqrt(g0./K.*(1+(K/km).^2));

    Lpm1 = exp( -5/4 * (kp./K).^2 );
    
    great_gama1 = exp(-(sqrt(K./kp)-1).^2./(2*sigma.^2));
    Jp1 = gama.^great_gama1;
    
    Bl1 = ap/2 * cp ./ cc1 .* Lpm1 .* Jp1 .* exp(-omega/sqrt(10)*(sqrt(K./kp)-1));
    %----------END 计算高频曲率谱----------%
    
    %----------Start 计算低频曲率谱----------%
    Bh1 = am/2*cm./cc1.*exp(-0.25*(K/km-1).^2);
    %----------END 计算低频曲率谱----------%
    
    %     Sk = K.^(-3).*(Bl + Bh);%全向谱Sk
    Sk1 = K.^(-3).*(Bl1+Lpm1.*Bh1);
    
    Fk1 = Sk1./K;% 计算二维谱使用的Fk
    
    [deltak1]=tony_spread(K,WindSpeed,omega);
    delta_phi=zeros(M+1,N+1);
    delta_phi(1:M/2,N/2+1:N+1)=pi;
    delta_phi(1:M/2,1:N/2)=pi;
    delta_phi(M/2+1:M+1,1:N/2)=2*pi;

    fine_transmit=delta_phi+atan((ky./kx));
    
    FD = Fk1./2/pi.*(1+deltak1.*cos(2.*((WindDirection-ObservingAngle)-fine_transmit)));
else
    error("Wrong ocean type input.");
end
end