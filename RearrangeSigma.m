%{
==========================================================================
Bacis Information
FileName:       RearrangeSigma.m
Author:         Shihao Zou
CreateTime:     2022.12
Version:        V1.0.0
Description:    This function is used for simulating the distance bunching
effect process.

==========================================================================
Input parameters
height: sea surface elevations.
dx: the size of spliting bin of sea surface along x axis.
dy: the size of spliting bin of sea surface along y axis.
Sigma0: backscattering coef of each grid of the sea surface.
Conf: Configuration file content

Output parameters:
Sigma0: Backscattering coef after distance bunching effect.
%}

function [Sigma0] = RearrangeSigma(height, Lx, Ly, dx, dy, Sigma0, Conf)
IncidentAngleIndex = Conf.ExpParameters.IncidentAngleIndex;
IncidentAngle = Conf.RadarParameters.IncidentAngles(IncidentAngleIndex) * pi/180;
RadarAltitude = Conf.RadarParameters.RadarAltitude;

[M, N] = size(height);

p = 0:M-1;
x = - Lx / 2 + p * dx;

radar_position_x = - RadarAltitude * tan(IncidentAngle);
radar_position_z = RadarAltitude;
looking_up_angle_tan = zeros(M, 1);
for i=1:M
    position_vec_x = radar_position_x - x(i);
    position_vec_z = radar_position_z - 0;
    looking_up_angle_tan(i) = abs(position_vec_x / position_vec_z);
end




end