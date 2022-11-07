function [deltak]=tony_spread(k,u10,omega)

CD10=1e-3*(0.8+0.065*u10);
us=sqrt(CD10)*u10;%%%%%%%%%%%%%%%%%%T.E���£�61��ʽ

cp=u10./omega;%%%%%%%%%%%%%%%%%%%��33��ʽ����һ��
g=9.81;
km=370.51;
c = sqrt(g ./ k); 
% c=sqrt(g./k.*(1+(k./km).^2));
cm=sqrt(2*g/km);
a0=log(2)/4;
ap=4;
am=0.13*us/cm;   %%%%%%%%%%%%%%%%��59��ʽ

deltak=tanh(a0+ap.*(c./cp).^2.5+am.*(cm./c).^2.5);%%%%%%%%%%%%%%%��57��ʽ

return

