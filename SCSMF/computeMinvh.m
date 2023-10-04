function [Minvh,const_hMh,c,b,stats] = computeMinvh(PD,FH)

Hk0=kron(PD.H,speye(PD.d^(PD.k0-1))); %make H to size of rho
[xindsInh,yindsInh]=Uinds(PD);

c = [1*PD.scale_rho ;zeros(PD.xinds(end,3) -1,1)]; 
b = PD.scale_sigma * [Hk0(:);zeros(PD.yinds(end,3)-PD.yinds(1,3),1)];
h=[c;b];
Minvh=nan(size(h));

[Minvh_x,Minvh_y,stats] = solveMinv(c,b,PD,FH,PD.MinvhTol,PD.MinvhMaxIter);

Minvh(xindsInh)= Minvh_x;
Minvh(yindsInh)= Minvh_y;

hTMinvh=h'*Minvh; 
const_hMh = 1/(1+hTMinvh);
end