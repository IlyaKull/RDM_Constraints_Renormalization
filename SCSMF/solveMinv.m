
function [z_x,z_y,stats] = solveMinv(w_x,w_y,PD,FH,tol,maxIter)
% see (28) in https://web.stanford.edu/~boyd/papers/pdf/scs_long.pdf

stats=struct('CGflag',[],'CGrelres',[],'CGiters',[],'T_CG',[]);
if nargin < 5
    tol=PD.CGtol;
end
if nargin < 6
   maxIter=PD.CGmaxIter;
end
    
ATw_y   = FH.ATransposed_funcHandle(full(w_y));
ticCG=tic;
[z_x,stats.CGflag,stats.CGrelres,stats.CGiters]     = pcg(FH.eyePlusATA_funcHandle,full(w_x - ATw_y),tol,maxIter);
stats.T_CG=toc(ticCG);
Az_x    = FH.A_funcHandle(z_x);
z_y     = w_y + Az_x;
