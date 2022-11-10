% Compute the following boundary layer paramaters: �߽����
% Cd10:drag coefficient, ����ϵ��
% ustar:friction velocity, Ħ���ٶ�
% zo:the roughness length,�ֲڶȳ���
% u195 and u125: the wind speed at 19.5 and 12.5 m height 
% Charnock parametrization is used

function [cd10, ustar, u195, u125, z0]=blpara(u10)

cd10 = 1E-3 * (0.8 + 0.065 * u10);    %initialization of cd10  ��T.E�����У�61��ʽ�е�Cd10N
ustar = sqrt(cd10) .* u10;       %%%%%%%%%%%%%T.E�����У�61��ʽ��ΪĦ������
const = 0.018;                    %Charnok constant
ck = 0.4;                         %Von Karmann constant
Va = 14E-6;                       % Air viscosity..................[m2.s]
g = 9.81;                         %Acceleration of gravity.........[m2/s]
it = 0;
ep = 1;
while ep > 0.01
    it = it + 1;
    z0 = 0.1 * Va ./ ustar + const * ustar .^ 2 / g;     %Roughness length  T.E�����У�64��ʽ
    cd10 = (ck ./ log(10 ./ z0) ) .^ 2;
    ustarn = sqrt(cd10) .* u10;
    ep = abs(ustar - ustarn) ./ ustar;
    ustar = ustarn;
end
u195 = (ustar / ck) .* log(19.5 ./ z0); 
u125 = (ustar / ck) .* log(12.5 ./ z0); 
