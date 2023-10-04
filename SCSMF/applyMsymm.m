function [Mxy] = applyMsymm(xy )
global PD
% applyMsymm takes  [x;y] and produces [x -A^T*y; -Ax-y];
[xindsInU,yindsInU]=Uinds();

Mxy=nan(size(xy)); %initialize
Mxy(xindsInU) =  xy(xindsInU) - applyAffineConstraintsTransposed(xy(yindsInU));
Mxy(yindsInU) = -xy(yindsInU) - applyAffineConstraints(xy(xindsInU));