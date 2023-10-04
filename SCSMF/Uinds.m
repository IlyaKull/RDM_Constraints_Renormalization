function [xindsInU,yindsInU,taoindInU]=Uinds(PD)
xindsInU=1:PD.xinds(end,3);
yindsInU=(1:PD.yinds(end,3)) + xindsInU(end);
taoindInU=yindsInU(end) +1;