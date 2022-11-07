%{
==========================================================================
Bacis Information
FileName:       AddingSpeckleNoise.m
Author:         Shihao Zou
CreateTime:     2018.04
Version:        V1.0.0
Description:    This function is used to add speckle noise into back soef
relativity. Note this method onyl holds true when the size of ocean
spliting bin equals to the resolution of given incident angle. For more
details about working principle please refer to the technique report 
"The Contruction Details of Forward and Inversion Software"

==========================================================================
Input parameters
Sigma0_Relativity: Relativity of back coef.
RadarGates： Radar gates number.


Output parameters:
Sigma0_noised: Relativity of back coef with speckle noise in space domain.
%}

function [ Sigma0_noised ] = AddingSpeckleNoise( Sigma0_Relativity, RadarGates )
k=RadarGates;
theta=Sigma0_Relativity/RadarGates;
Sigma0_noised=gamrnd(k,theta);
end