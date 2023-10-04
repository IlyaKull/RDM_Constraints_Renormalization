function [eyePlusQu] = eyePlusQ(u,PD,FH)
% apply eye+Q (see https://web.stanford.edu/~boyd/papers/pdf/scs_long.pdf 4.1)
% eye+Q = [M,h;-h^T,1] ie (eye+Q)[xy,tao] = [M*[xy]+h*tao ; -h^T*xy + tao]
 

[xindsInU,yindsInU,taoindInU]=Uinds(PD);
xyindsInU=[xindsInU,yindsInU];

 
h=[PD.c;PD.b];


eyePlusQu=nan(size(u)); %initialize
eyePlusQu(xyindsInU) =  applyM(u(xyindsInU),PD,FH) + h*u(taoindInU);
eyePlusQu(taoindInU) =  -h'*u(xyindsInU) + u(taoindInU);