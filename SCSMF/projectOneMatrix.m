function [P,asymm]= projectOneMatrix(S)
symmS=0.5*(S+S');
asymm = norm(S-S');

[V,D] = eig(symmS,'vector');
 
D(D<0)=0;
P=real(V*diag(D)*V');