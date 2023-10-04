function [u,v,stats]=iterateSCS(u,v,iter,PD,FH)
% see https://web.stanford.edu/~boyd/papers/pdf/scs_long.pdf 3.2.3 and 3.3
% q is the over-relaxation parameter PD.q
 
ticAff=tic;
[utilde,stats]= projectToAffine(u+v,iter,PD,FH);
stats.T_Aff=toc(ticAff);
qcomb   = (PD.q*utilde + (1-PD.q)*u );

ticProj=tic;
[u , stats.avAsymm] = projectToCone(qcomb - v,PD);
stats.T_Proj=toc(ticProj);

v = v - qcomb + u; 