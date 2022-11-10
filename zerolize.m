function [amn, F] = zerolize(amn, F, fine_D, fine_thistype, kx_odds, ky_odds )

%-------------1.给定单边波高谱朝向，并确定其分割线方程
%-------------划分线集合为{ky=tan(phi)*kx,kx=0},其中phi取值范围为（pi/2,2*pi+pi/2）
theta = mod(fine_thistype - fine_D, 2 * pi);
phi = theta + pi / 2;%phi是分割线的与kx正方向的夹角，与波高谱朝向角theta相差90°

%-------------2.根据波高谱朝向，修改kx_odds,ky_odds,k_combined矩阵
%-------------引入判据矩阵atr，代表分割线ky值
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
%---------------------------END 对称的波高谱置零得到单边海波高面---------------------------%
end