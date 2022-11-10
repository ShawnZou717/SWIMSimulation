function [amn, F] = zerolize(amn, F, fine_D, fine_thistype, kx_odds, ky_odds )

%-------------1.�������߲����׳��򣬲�ȷ����ָ��߷���
%-------------�����߼���Ϊ{ky=tan(phi)*kx,kx=0},����phiȡֵ��ΧΪ��pi/2,2*pi+pi/2��
theta = mod(fine_thistype - fine_D, 2 * pi);
phi = theta + pi / 2;%phi�Ƿָ��ߵ���kx������ļнǣ��벨���׳����theta���90��

%-------------2.���ݲ����׳����޸�kx_odds,ky_odds,k_combined����
%-------------�����оݾ���atr������ָ���kyֵ
if theta == 0 || theta == 2 * pi
    amn(kx_odds <= 0) = 0;
    F(kx_odds <= 0) = 0;
elseif theta == pi
    amn(kx_odds >= 0) = 0;
    F(kx_odds >= 0) = 0;
elseif 0 < theta && theta < pi
    atr = tan(phi) .* kx_odds;
    id_atr = find(atr >= ky_odds);
    amn(id_atr) = 0;
    F(id_atr) = 0;
elseif pi < theta && theta < 2 * pi
    atr = tan(phi) .* kx_odds;
    id_atr = find(atr <= ky_odds);
    amn(id_atr) =0;
    F(id_atr) = 0; 
end
%---------------------------END �ԳƵĲ���������õ����ߺ�������---------------------------%
end