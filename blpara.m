% Compute the following boundary layer paramaters: 边界参数
% Cd10:drag coefficient, 阻力系数
% ustar:friction velocity, 摩擦速度
% zo:the roughness length,粗糙度长度
% u195 and u125: the wind speed at 19.5 and 12.5 m height 
% Charnock parametrization is used

function [cd10, ustar, u195, u125, z0]=blpara(u10)

cd10 = 1E-3 * (0.8 + 0.065 * u10);    %initialization of cd10  即T.E文章中（61）式中的Cd10N
ustar = sqrt(cd10) .* u10;       %%%%%%%%%%%%%T.E文章中（61）式，为摩擦因子
const = 0.018;                    %Charnok constant
ck = 0.4;                         %Von Karmann constant
Va = 14E-6;                       % Air viscosity..................[m2.s]
g = 9.81;                         %Acceleration of gravity.........[m2/s]
it = 0;
ep = 1;
while ep > 0.01
    it = it + 1;
    z0 = 0.1 * Va ./ ustar + const * ustar .^ 2 / g;     %Roughness length  T.E文章中（64）式
    cd10 = (ck ./ log(10 ./ z0) ) .^ 2;
    ustarn = sqrt(cd10) .* u10;
    ep = abs(ustar - ustarn) ./ ustar;
    ustar = ustarn;
end
u195 = (ustar / ck) .* log(19.5 ./ z0); 
u125 = (ustar / ck) .* log(12.5 ./ z0); 
