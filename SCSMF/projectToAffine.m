function [u,stats] = projectToAffine(w,iter,PD,FH)
% solves (I+Q)u=w
% see https://web.stanford.edu/~boyd/papers/pdf/scs_long.pdf (4.1)
% PD is the problem data structure. this function adds the Minvh field and
% others
% if Minv*h is not cached compute it  
 
[xindsInU,yindsInU,taoindInU]=Uinds(PD);

% h= [c;b];
w_tao=w(taoindInU);
w_xIN=w(xindsInU) -w_tao*(PD.c);
w_yIN=w(yindsInU) -w_tao*(PD.b);

[out_x,out_y,stats] = solveMinv(w_xIN ,w_yIN,PD,FH, PD.CGtolFunc(iter));
out=[out_x;out_y];
u_xy = out - PD.const_hMh * PD.Minvh * ([PD.c;PD.b]' * out);
% ad
u = [u_xy; w_tao + (PD.c)'*u_xy(xindsInU) + (PD.b)'*u_xy(yindsInU)];

% stats.testSol=norm(w-eyePlusQ(u));
end


