%{
==========================================================================
Bacis Information
FileName:       GeneratingSeaSurface.m
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
t:  Current time. (t second)
ObservingAngle: Angle between scanning direction and the flight direction.

Output parameters:
height: sea surface elevation.
slope_x: slope of sea surface along x axis.
slope_y: slope of sea surface along y axis.
%}

function [ height, slope_x, slope_y ] = GeneratingSeaSurface( Lx, Ly, dx, dy, t, ObservingAngle )
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


% index matrix in Cartesian space domain
p=m;
x = -Lx/2+p*dx;
q=n;
y = -Ly/2+q*dx;

% generating initial phase matrix
init_phase = GeneratingInitPhase(M,N);
% Angular frequency accroding to dispersion relation
g0 = 9.81;
wn=(K*g0).^0.5;

Conf = jsondecode(fileread("ParaConfig.json"));
OceanWaveType = Conf.ExpParameters.OceanWaveType;
OceanParameters = Conf.ExpParameters.OceanParameters;
L = length(OceanWaveType);

% initial wave spectrum
F = 0;
% initial amplitude
amn = 0;

% Generating zerolized wave spectrum. For more details please refer to the 
% Patent of Invention "A Monte-Carlo-based method to generate one-way 
% propagating ocean waves"
for H=1:L
    ocean_type = OceanWaveType(H);
    ocean_paras = OceanParameters(H);
    ocean_paras = ocean_paras{1,1};
    [F0] = SpectrumConstruction(Lx, Ly, dx, dy, ObservingAngle, ocean_type, ocean_paras);

    amn0=2*sqrt(F0.*dkx.*dky);
    amn0(isnan(amn0))=0;

    [ amn0,F0 ] = zerolize( amn0, F0, ObservingAngle, ocean_paras.PeakWaveDirection, kx, ky );
    amn = amn+amn0;
    F = F+F0;
end

% Details about the principle of sign_matrix can be found in the patent of
% invention "A Monte-Carlo-based method to generate one-way propagating 
% ocean waves"
sign_matrix = exp(1i*pi*(m+n));
middle_term = init_phase.*amn.*sign_matrix.*exp(1i.*(-wn*t));
sign_matrix = sign_matrix(2:end,2:end);
middle_term = middle_term(2:end,2:end);

middle_sx = middle_term.*kx(2:end,2:end)*1i;
middle_sy = middle_term.*ky(2:end,2:end)*1i;

height=M*N*ifft2(middle_term).*sign_matrix;
slope_x=M*N*ifft2(middle_sx).*sign_matrix;
slope_y=M*N*ifft2(middle_sy).*sign_matrix;

height=real(height);
% Surely the slope coule be calculated just by differential height, but
% that would introduce estimation error.
slope_x=real(slope_x);
slope_y=real(slope_y);

end