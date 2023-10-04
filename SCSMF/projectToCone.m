function [Pu,avAsymm] = projectToCone(u,PD)
% u = [ x y tao]^T
% the projection is to ( R^n x PSD x R+ ) ie
% x maps to itself, y gets projected to the positive semidefinite cone
% elementwise (each matrix composing y gets projected); and tao is projected to
% R+

Pu=u; % initialize

% x goes to x (projection to R^n)

% project y
[~,yindsInU]=Uinds(PD);
[Pu(yindsInU),avAsymm]=projectToPSDcone(u(yindsInU),PD);

% project tao to R_+
if Pu(end)<0
    Pu(end)=0;
end
