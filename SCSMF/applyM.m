function [Mxy] = applyM(xy,PD,FH)
% M= [I,A^T; -A,I]
% applyMsymm takes  [x;y] and produces [x +A^T*y; -Ax+y];
 
[xindsInU,yindsInU]=Uinds(PD);

Mxy=nan(size(xy)); %initialize
Mxy(xindsInU) = xy(xindsInU) + FH.ATransposed_funcHandle(xy(yindsInU));
Mxy(yindsInU) = xy(yindsInU) - FH.A_funcHandle(xy(xindsInU));