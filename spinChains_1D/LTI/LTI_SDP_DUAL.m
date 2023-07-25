function [ ElocTIRig,ElocTI] = LTI_SDP_DUAL(N,H,sdpOps)
%solve the locally translation invariant problem for a state of N spins

if isstruct(H)
    suppH=H.suppH; % the number of sites H acts on 
    d = H.d; % the dimention of each site
    H=H.Hamiltonian; % the Hamiltonian itself
else % H is a d^2 x d^2 matrix
    suppH=2;
    d=sqrt(size(H,1));
end
 
 
alpha=sdpvar(d^(N-1)); %first state is on k0+1 sites. Alpha is the local transl. inv. dual variable

sdpvar epsil

expr = kron(H,speye(d^(N-suppH))) ...
    + dontUseKron(alpha,d^(N-1),d,'IX') ...
    - dontUseKron(alpha,d^(N-1),d,'XI') ...
    - epsil*speye(d^(N)); 

constraints = [ expr >= 0 ];
 
objectiveFunc = - epsil; %minus sign for maximization

 optimize(constraints,objectiveFunc,sdpOps);
% %%
ElocTI=-value(objectiveFunc);
 
expr = value(expr);
minEig = min(real(eig( expr )));
 
ElocTIRig=value(epsil) + minEig;
 