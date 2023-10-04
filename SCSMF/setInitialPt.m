function [u,v]=setInitialPt(PD)

[~,~,taoindInU]=Uinds(PD);
u=zeros(taoindInU,1);
u(taoindInU)=1;
v=u;
