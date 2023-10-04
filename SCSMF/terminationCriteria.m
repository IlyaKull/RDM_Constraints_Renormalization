function [energy,primalRes,dualRes,gapRes,tolReached,x,y,s,vMinusQu]=terminationCriteria(u,v,PD,FH)
% see https://web.stanford.edu/~boyd/papers/pdf/scs_long.pdf  (3.5)
% u=[x,y,tao]'  ; v = [r,s,kappa]'
[xindsInU,yindsInU]=Uinds(PD);
c_scaled=PD.c/PD.scale_rho;
b_scaled=PD.b/PD.scale_sigma;


x=u(xindsInU)/u(end)/PD.scale_sigma;
y=u(yindsInU)/u(end)/PD.scale_rho;
s=v(yindsInU)/u(end)/PD.scale_sigma;
energy=-(c_scaled'*x);

% primal residuals p=Ax+s-b, !! x=u_x/u_tao , s=v_s/u_tao
primalRes = norm(FH.A_funcHandle(x) + s - b_scaled);

% dual residuals d=A^T*y+c !! y=u_y/u_tao
dualRes = norm(FH.ATransposed_funcHandle(y) + c_scaled);

% duality gap tolerance 
gapRes = -energy + b_scaled'*y;

if nargout >8
    vMinusQu=norm(eyePlusQ(u,PD)-u-v);
end

tolReached= all([primalRes, dualRes, abs(gapRes)] ./PD.residTol <=  [1,1,1] ) ;
           
